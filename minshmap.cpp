// minSHmap (C++) - a line-for-line port of minshmap.py: same algorithm, byte-identical
// output. The only extra is `-j` (parallel reads). Keep it in lockstep with minshmap.py.
// sketch: FracMinHash rolling ntHash | index: hash->[(seg,pos,strand)] | map: rank read
// k-mers rarest-first, scatter rarest into overlapping windows, score containment, prune
// via sh = 1-(used-matches)/m. Build: g++ -O3 -std=c++17 -march=native -pthread minshmap.cpp
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>
using u64 = uint64_t;
using std::string;
using std::vector;
static const u64 MASK64 = ~u64(0);
static u64 LUT_FW[256], LUT_RC[256];  // a k-mer and its rev-comp share a hash
static void init_lut() {
    LUT_FW['A'] = 0x3C8BFBB395C60474ULL; LUT_FW['C'] = 0x3193C18562A02B4CULL;
    LUT_FW['G'] = 0x20323ED082572324ULL; LUT_FW['T'] = 0x295549F54BE24456ULL;
    LUT_RC['A'] = LUT_FW['T']; LUT_RC['C'] = LUT_FW['G'];
    LUT_RC['G'] = LUT_FW['C']; LUT_RC['T'] = LUT_FW['A'];
}
static inline u64 rotl(u64 x, int r) { r &= 63; return r ? (x << r) | (x >> (64 - r)) : x; }
static inline u64 rotr(u64 x, int r) { r &= 63; return r ? (x >> r) | (x << (64 - r)) : x; }
struct SkEntry { int pos; u64 h; int strand; };  // strand: 0 forward, 1 rev-comp

// Kept (pos, hash, strand) k-mers; mirrors minshmap.py sketch(). O(len).
static vector<SkEntry> sketch(const string &seq, int k, double hfrac) {
    vector<SkEntry> out;
    int n = (int)seq.size();
    if (n < k) return out;
    u64 thr = hfrac >= 1.0 ? MASK64 : (u64)(hfrac * 18446744073709551616.0), h_fw = 0, h_rc = 0;
    for (int t = 0; t < k; ++t) {                      // hash of the first window
        h_fw ^= rotl(LUT_FW[(unsigned char)seq[t]], k - 1 - t);
        h_rc ^= rotl(LUT_RC[(unsigned char)seq[t]], t);
    }
    for (int i = 0; i + k <= n; ++i) {
        u64 h = h_fw <= h_rc ? h_fw : h_rc;            // canonical: min over both strands
        if (h <= thr) out.push_back({i, h, h_fw <= h_rc ? 0 : 1});
        if (i + k < n) {                               // roll the window one base
            unsigned char o = seq[i], in = seq[i + k];
            h_fw = rotl(h_fw, 1) ^ rotl(LUT_FW[o], k) ^ LUT_FW[in];
            h_rc = rotr(h_rc, 1) ^ rotr(LUT_RC[o], 1) ^ rotl(LUT_RC[in], k - 1);
        }
    }
    return out;
}
struct Hit { int sid; int pos; int strand; };
struct Segment { string name; int length; };
struct Index { std::unordered_map<u64, vector<Hit>> h2hits; vector<Segment> segments; };

static Index build_index(const vector<std::pair<string, string>> &refs, int k, double hfrac) {
    Index idx;                                         // hash -> [(segm_id, pos, strand), ...]
    for (size_t sid = 0; sid < refs.size(); ++sid) {
        idx.segments.push_back({refs[sid].first, (int)refs[sid].second.size()});
        for (const auto &e : sketch(refs[sid].second, k, hfrac))
            idx.h2hits[e.h].push_back({(int)sid, e.pos, e.strand});
    }
    return idx;
}
struct Mapping { int sid; int t_start; int t_end; double score; int codir; bool ok = false; };

// Best reference window for the read, or ok=false. Mirrors minshmap.py map_read().
static Mapping map_read(const string &seq, const Index &idx, int k, double hfrac, double theta) {
    Mapping best;
    auto sk = sketch(seq, k, hfrac);
    int m = (int)sk.size();                            // informative k-mers in the read
    if (m == 0) return best;
    int W = std::max((int)seq.size(), 1);              // candidate windows are read-length-wide
    struct Seed { const vector<Hit> *hits; int n_hits; int rstrand; };
    std::unordered_set<u64> seen;
    vector<Seed> seeds;
    for (const auto &e : sk) {                         // one seed per distinct read k-mer
        if (!seen.insert(e.h).second) continue;
        auto it = idx.h2hits.find(e.h);
        const vector<Hit> *hits = it == idx.h2hits.end() ? nullptr : &it->second;
        seeds.push_back({hits, hits ? (int)hits->size() : 0, e.strand});
    }
    std::stable_sort(seeds.begin(), seeds.end(),       // rarest first (stable: ties keep order)
                     [](const Seed &a, const Seed &b) { return a.n_hits < b.n_hits; });
    int S = (int)((1.0 - theta) * m) + 1;              // this many rare seeds are enough
    std::unordered_set<long long> cand;                // window keys: (sid<<32)|bucket
    for (int i = 0; i < S && i < (int)seeds.size(); ++i)
        if (seeds[i].hits)
            for (const Hit &hit : *seeds[i].hits) {    // overlapping buckets b and b-1
                int b = hit.pos / W;
                cand.insert(((long long)hit.sid << 32) | (unsigned)b);
                if (b > 0) cand.insert(((long long)hit.sid << 32) | (unsigned)(b - 1));
            }
    vector<long long> order(cand.begin(), cand.end());
    std::sort(order.begin(), order.end());             // deterministic, matches Python sorted(cand)
    for (long long key : order) {                      // score each window, keep the best
        int sid = (int)(key >> 32), b = (int)(unsigned)(key & 0xffffffffLL);
        int lo = b * W, hi = (b + 2) * W, used = 0, matches = 0, codir = 0, r_min = -1, r_max = -1;
        bool pruned = false;
        for (const Seed &s : seeds) {                  // add seeds rarest-first, prune early
            ++used;
            if (s.hits)
                for (const Hit &hit : *s.hits)
                    if (hit.sid == sid && lo <= hit.pos && hit.pos < hi) {  // count once per k-mer
                        ++matches;
                        codir += (hit.strand == s.rstrand) ? 1 : -1;
                        r_min = r_min < 0 ? hit.pos : std::min(r_min, hit.pos);
                        r_max = std::max(r_max, hit.pos);
                        break;
                    }
            if (1.0 - double(used - matches) / m < theta) { pruned = true; break; }  // SEED HEURISTIC
        }
        double score = double(matches) / m;            // containment = fraction of read k-mers hit
        if (!pruned && score >= theta && (!best.ok || score > best.score))
            best = {sid, r_min, r_max + k, score, codir, true};
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
    bool have = false;
    while (std::getline(f, line)) {
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
                      const Index &idx, int k, double hfrac, double theta, vector<string> &out) {
    for (size_t i = lo; i < hi; ++i) {
        Mapping best = map_read(reads[i].second, idx, k, hfrac, theta);
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
    init_lut();
    int k = 15, threads = 1;
    double hfrac = 0.05, theta = 0.9;
    vector<string> pos;
    for (int i = 1; i < argc; ++i) {
        string a = argv[i];
        auto next = [&]() { return string(argv[++i]); };
        if (a == "-k") k = std::stoi(next());
        else if (a == "-r" || a == "--hfrac") hfrac = std::stod(next());
        else if (a == "-t" || a == "--theta") theta = std::stod(next());
        else if (a == "-j" || a == "--threads") threads = std::stoi(next());
        else if (!a.empty() && a[0] != '-') pos.push_back(a);
        else { std::cerr << "Unknown option: " << a << "\n"; return 1; }
    }
    if (pos.size() < 2) { std::cerr << "Usage: minshmap ref.fa reads.fa [-k][-r][-t][-j]\n"; return 1; }
    Index idx = build_index(read_fasta(pos[0]), k, hfrac);
    auto reads = read_fasta(pos[1]);
    vector<string> out(reads.size());
    parallel_for(reads.size(), threads, [&](size_t lo, size_t hi) { map_chunk(reads, lo, hi, idx, k, hfrac, theta, out); });
    for (const string &line : out) std::cout << line;   // original read order preserved
    return 0;
}
