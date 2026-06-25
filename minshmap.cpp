// minSHmap (C++) - sketch-based read mapper; same algorithm as the pedagogical minshmap.py,
// optimized: flat CSR hit layout, ankerl::unordered_dense maps, reserve(), `-j` parallel reads.
// Byte-identical output to minshmap.py (same hash, map order = sorted(cand)). Pipeline: (w,k)-minimizers
// (2-bit rolling code + SplitMix64) -> hash->[(seg,pos,strand)] index -> map: rank read minimizers
// rarest-first, scatter into overlapping windows, score containment, prune sh=1-(used-matches)/m.
// Build: g++ -O3 -std=c++17 -march=native -pthread -I ../shmap/ext/unordered_dense/include -o minshmap minshmap.cpp
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <deque>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <thread>
#include <vector>
#include "ankerl/unordered_dense.h"
using u64 = uint64_t;
using std::string;
using std::vector;
static const u64 MASK64 = ~u64(0);
static inline u64 mix(u64 x) {                          // SplitMix64 finalizer: packed code -> 64-bit hash
    x += 0x9E3779B97F4A7C15ULL; x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL; return x ^ (x >> 31);
}
static inline int b2(unsigned char c) { return c == 'C' ? 1 : c == 'G' ? 2 : c == 'T' ? 3 : 0; }  // complement(b)=3-b
struct SkEntry { int pos; u64 h; int strand; };  // strand: 0 forward, 1 rev-comp

// (w,k)-minimizers; mirrors minshmap.py minimizers(). Each window of w consecutive k-mers
// contributes its smallest-hash k-mer (leftmost on ties), emitted once in pos order. k<=32.
static vector<SkEntry> sketch(const string &seq, int k, int w) {
    vector<SkEntry> out;
    int n = (int)seq.size();
    if (n < k) return out;
    u64 mask = k >= 32 ? MASK64 : ((u64(1) << (2 * k)) - 1), fw = 0, rc = 0;    // window mask + fw/rc 2-bit codes
    int sh = 2 * (k - 1);
    auto roll = [&](int b) { fw = ((fw << 2) | b) & mask; rc = (rc >> 2) | ((u64)(3 - b) << sh); };
    for (int t = 0; t < k; ++t) roll(b2(seq[t]));       // pack the first window
    vector<SkEntry> km;                                 // (pos, hash, strand) for every k-mer
    for (int i = 0;; ++i) {
        u64 c = fw <= rc ? fw : rc;                     // canonical k-mer (min over both strands)
        km.push_back({i, mix(c), fw <= rc ? 0 : 1});
        if (i + k >= n) break;
        roll(b2(seq[i + k]));                            // roll the window one base
    }
    int K = (int)km.size(), ww = std::min(w, K), last = -1;
    std::deque<int> dq;                                 // indices; hashes increasing, front = window min
    for (int t = 0; t < K; ++t) {                       // sliding-window minimum, O(K)
        while (!dq.empty() && km[dq.back()].h > km[t].h) dq.pop_back();  // strict -> leftmost on ties
        dq.push_back(t);
        if (dq.front() <= t - ww) dq.pop_front();        // evict indices left of the window
        if (t >= ww - 1 && dq.front() != last) out.push_back(km[last = dq.front()]);
    }
    return out;
}
struct Hit { int sid; int pos; int strand; };
struct Segment { string name; int length; };
// CSR layout: hits for a hash are contiguous in `hits`; `id` maps hash -> dense slot s,
// spanning hits[off[s]..off[s+1]). One allocation (not a vector per hash) -> fewer cache
// misses; hits in a slot stay sorted by (sid, pos) so map_read can binary-search them.
struct Index {
    ankerl::unordered_dense::map<u64, uint32_t> id;
    vector<uint32_t> off;
    vector<Hit> hits;
    vector<Segment> segments;
    const Hit *range(u64 h, int &n) const {            // hits for hash h (n = count)
        auto it = id.find(h);
        if (it == id.end()) { n = 0; return nullptr; }
        uint32_t s = it->second;
        n = (int)(off[s + 1] - off[s]);
        return hits.data() + off[s];
    }
};

static Index build_index(const vector<std::pair<string, string>> &refs, int k, int w) {
    Index idx;
    vector<std::pair<u64, Hit>> raw;                   // (hash, hit) in reference order
    for (size_t sid = 0; sid < refs.size(); ++sid) {
        idx.segments.push_back({refs[sid].first, (int)refs[sid].second.size()});
        for (const auto &e : sketch(refs[sid].second, k, w))
            raw.push_back({e.h, {(int)sid, e.pos, e.strand}});
    }
    idx.id.reserve(raw.size());
    vector<uint32_t> slot(raw.size()), cnt;            // dense slot per entry; hits per slot (first-seen order)
    for (size_t i = 0; i < raw.size(); ++i) {          // single hashmap lookup/entry, remember the slot
        uint32_t s = idx.id.try_emplace(raw[i].first, (uint32_t)cnt.size()).first->second;
        if (s == (uint32_t)cnt.size()) cnt.push_back(0);
        ++cnt[s]; slot[i] = s;
    }
    idx.off.assign(cnt.size() + 1, 0);
    for (size_t s = 0; s < cnt.size(); ++s) idx.off[s + 1] = idx.off[s] + cnt[s];
    idx.hits.resize(raw.size());
    vector<uint32_t> cur(idx.off.begin(), idx.off.end() - 1);  // running write head per slot
    for (size_t i = 0; i < raw.size(); ++i) idx.hits[cur[slot[i]]++] = raw[i].second;  // stable placement
    return idx;
}
struct Mapping { int sid; int t_start; int t_end; double score; int codir; bool ok = false; };
struct Seed { const Hit *hits; int n_hits; int rstrand; };  // one distinct read minimizer

// Upper bound on the containment a window can still reach, in [0,1] (compare to theta).
static inline double seed_heuristic(int used, int matches, int m) { return 1.0 - double(used - matches) / m; }

// One seed per distinct read minimizer (its ref hits + read strand), sorted rarest-first.
static vector<Seed> gather_seeds(const vector<SkEntry> &sk, const Index &idx) {
    ankerl::unordered_dense::set<u64> seen;
    vector<Seed> seeds;
    seen.reserve(sk.size()); seeds.reserve(sk.size());
    for (const auto &e : sk) {
        if (!seen.insert(e.h).second) continue;
        int nh; const Hit *hp = idx.range(e.h, nh);
        seeds.push_back({hp, nh, e.strand});
    }
    std::stable_sort(seeds.begin(), seeds.end(),       // rarest first (stable: ties keep order)
                     [](const Seed &a, const Seed &b) { return a.n_hits < b.n_hits; });
    return seeds;
}

// Scatter the rarest seeds' hits into overlapping buckets (b and b-1); window keys by votes desc.
static vector<long long> candidate_windows(const vector<Seed> &seeds, int m, double theta, int W) {
    int S = (int)((1.0 - theta) * m) + 1;              // this many rare seeds are enough
    ankerl::unordered_dense::map<long long, int> cand; // key (sid<<32)|bucket -> votes
    for (int i = 0; i < S && i < (int)seeds.size(); ++i)
        for (int j = 0; j < seeds[i].n_hits; ++j) {
            const Hit &hit = seeds[i].hits[j];
            int b = hit.pos / W;
            ++cand[((long long)hit.sid << 32) | (unsigned)b];
            if (b > 0) ++cand[((long long)hit.sid << 32) | (unsigned)(b - 1)];
        }
    vector<std::pair<int, long long>> order;
    order.reserve(cand.size());
    for (const auto &kv : cand) order.push_back({kv.second, kv.first});
    std::sort(order.begin(), order.end(), [](const auto &a, const auto &b) {
        return a.first != b.first ? a.first > b.first : a.second < b.second;  // votes desc, key asc
    });
    vector<long long> keys;
    keys.reserve(order.size());
    for (const auto &o : order) keys.push_back(o.second);
    return keys;
}
struct WinScore { int matches, codir, r_min, r_max; bool pruned; };

// Add seeds rarest-first to window [lo,hi); stop once seed_heuristic proves it can't reach
// target. Each seed contributes its smallest-pos hit (binary-searched in the sorted slot).
static WinScore score_window(const vector<Seed> &seeds, int sid, int lo, int hi, int m, double target) {
    int used = 0, matches = 0, codir = 0, r_min = -1, r_max = -1;
    for (const Seed &s : seeds) {
        ++used;
        int a = 0, z = s.n_hits;                        // lower_bound (sid,lo): smallest-pos match, O(log n)
        while (a < z) {
            int mid = (a + z) >> 1;
            const Hit &h = s.hits[mid];
            if (h.sid < sid || (h.sid == sid && h.pos < lo)) a = mid + 1;
            else z = mid;
        }
        if (a < s.n_hits && s.hits[a].sid == sid && s.hits[a].pos < hi) {
            const Hit &hit = s.hits[a]; ++matches;
            codir += (hit.strand == s.rstrand) ? 1 : -1;
            r_min = r_min < 0 ? hit.pos : std::min(r_min, hit.pos);
            r_max = std::max(r_max, hit.pos);
        }
        if (seed_heuristic(used, matches, m) < target) return {matches, codir, r_min, r_max, true};
    }
    return {matches, codir, r_min, r_max, false};
}

// Best reference window for the read, or ok=false. Mirrors minshmap.py map_read().
static Mapping map_read(const string &seq, const Index &idx, int k, int w, double theta) {
    Mapping best;
    auto sk = sketch(seq, k, w);
    int m = (int)sk.size();                            // informative minimizers in the read
    if (m == 0) return best;
    int W = std::max((int)seq.size(), 1);              // candidate windows are read-length-wide
    vector<Seed> seeds = gather_seeds(sk, idx);
    vector<long long> windows = candidate_windows(seeds, m, theta, W);
    long long best_key = 0;                            // tie-break on key -> order-independent result
    for (long long key : windows) {                    // score promising windows first
        int sid = (int)(key >> 32), b = (int)(unsigned)(key & 0xffffffffLL);
        double target = best.ok ? std::max(theta, best.score) : theta;  // can't beat best -> prune harder
        WinScore ws = score_window(seeds, sid, b * W, (b + 2) * W, m, target);
        if (ws.pruned) continue;
        double score = double(ws.matches) / m;         // containment = fraction of read minimizers hit
        if (score >= theta && (!best.ok || score > best.score || (score == best.score && key < best_key))) {
            best = {sid, ws.r_min, ws.r_max + k, score, ws.codir, true};
            best_key = key;
        }
    }
    return best;
}
static string first_token(const string &s) { size_t p = s.find_first_of(" \t"); return p == string::npos ? s : s.substr(0, p); }

// (name, sequence) per FASTA record; name is the first header token.
static vector<std::pair<string, string>> read_fasta(const string &path) {
    vector<std::pair<string, string>> recs;
    std::ifstream f(path);
    if (!f) { std::cerr << "Cannot open " << path << "\n"; std::exit(1); }
    string line, name, seq;
    bool have = false;    while (std::getline(f, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (have) recs.push_back({name, seq});
            name = first_token(line.substr(1)); seq.clear(); have = true;
        } else { for (char &c : line) c = (char)std::toupper((unsigned char)c); seq += line; }
    }
    if (have) recs.push_back({name, seq});
    return recs;
}

// Map reads[lo,hi) into out[lo,hi); threads write disjoint slices, no locks -> byte-identical.
static void map_chunk(const vector<std::pair<string, string>> &reads, size_t lo, size_t hi,
                      const Index &idx, int k, int w, double theta, vector<string> &out) {
    for (size_t i = lo; i < hi; ++i) {
        Mapping best = map_read(reads[i].second, idx, k, w, theta);
        if (!best.ok) continue;
        const Segment &seg = idx.segments[best.sid];
        int qlen = (int)reads[i].second.size();
        int nmatch = (int)std::rint(best.score * qlen);   // round-half-to-even, like Python round()
        std::ostringstream os;
        os << reads[i].first << '\t' << qlen << '\t' << 0 << '\t' << qlen << '\t'
           << (best.codir >= 0 ? '+' : '-') << '\t' << seg.name << '\t' << seg.length << '\t'
           << best.t_start << '\t' << best.t_end << '\t' << nmatch << '\t'
           << (best.t_end - best.t_start) << '\t' << 60 << '\n';   // mapq constant 60
        out[i] = os.str();
    }
}

// Run work(lo, hi) over [0, n) in `threads` contiguous stripes; <=1 runs inline.
template <class F>
static void parallel_for(size_t n, int threads, F work) {
    if (threads <= 1 || n < 2) { work(0, n); return; }
    threads = (int)std::min<size_t>(threads, n);
    vector<std::thread> pool;
    size_t chunk = (n + threads - 1) / threads;
    for (int t = 0; t < threads; ++t) {
        size_t lo = std::min(n, (size_t)t * chunk), hi = std::min(n, lo + chunk);
        if (lo < hi) pool.emplace_back(work, lo, hi);
    }
    for (auto &th : pool) th.join();
}

int main(int argc, char **argv) {
    int k = 15, w = 10, threads = 1;
    double theta = 0.9;
    vector<string> pos;
    for (int i = 1; i < argc; ++i) {
        string a = argv[i];
        auto next = [&]() { return string(argv[++i]); };
        if (a == "-k") k = std::stoi(next());
        else if (a == "-w" || a == "--window") w = std::stoi(next());
        else if (a == "-t" || a == "--theta") theta = std::stod(next());
        else if (a == "-j" || a == "--threads") threads = std::stoi(next());
        else if (!a.empty() && a[0] != '-') pos.push_back(a);
        else { std::cerr << "Unknown option: " << a << "\n"; return 1; }
    }
    if (pos.size() < 2) { std::cerr << "Usage: minshmap ref.fa reads.fa [-k][-w][-t][-j]\n"; return 1; }
    Index idx = build_index(read_fasta(pos[0]), k, w);
    auto reads = read_fasta(pos[1]);
    vector<string> out(reads.size());
    parallel_for(reads.size(), threads, [&](size_t lo, size_t hi) { map_chunk(reads, lo, hi, idx, k, w, theta, out); });
    for (const string &line : out) std::cout << line;   // original read order preserved
    return 0;
}
