// minSHmap (C++) - optimized port of minshmap.py (CSR index, binary-searched blocks, `-j` reads).
// Minimizers come from the SAME library as the Python tool: rust-seq/minimizer-iter via the
// minimizer_ext C ABI (so py and cpp stay byte-identical, no copy-pasted sketch). w must be ODD.
// Reports a phi-FREE mapq (Def. 5'/6'): an alternative is any scored block whose reference
// interval is DISJOINT from the best; mapq=60 iff every alternative is weaker by > delta.
// Build: cargo build --release --no-default-features (in minimizer_ext/), then
//   g++ -O3 -std=c++17 -march=native -pthread -I ../shmap/ext/unordered_dense/include -o minshmap minshmap.cpp \
//       -L minimizer_ext/target/release -lminimizer_ext -lws2_32 -luserenv -lbcrypt -lntdll
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <thread>
#include <vector>
#include "ankerl/unordered_dense.h"
#include "bench.hpp"                        // optional timing/memory (only when MINSHMAP_BENCH is set)
using u64 = uint64_t; using std::string; using std::vector;
namespace adh = ankerl::unordered_dense;
extern "C" { struct Mz { u64 pos, val; uint8_t strand; }; Mz *mz_compute(const char *seq, size_t len, size_t k, size_t w, size_t *n); void mz_free(Mz *p, size_t n); }
struct SkEntry { int pos; u64 h; int strand; };
struct Hit { uint32_t sid, pv; int pos() const { return (int)(pv >> 1); } int strand() const { return (int)(pv & 1u); } };  // 8B: pv=(pos<<1)|strand
static inline Hit make_hit(int sid, int pos, int strand) { return { (uint32_t)sid, (uint32_t)(((uint32_t)pos << 1) | (uint32_t)(strand & 1)) }; }
struct Segment { string name; int length; };
// (w,k)-minimizers straight from minimizer-iter (rust-seq), the exact same lib minshmap.py uses. w must be ODD.
static vector<SkEntry> sketch(const string &seq, int k, int w) {
    vector<SkEntry> out; size_t n = 0;
    Mz *mz = mz_compute(seq.data(), seq.size(), (size_t)k, (size_t)w, &n);
    if (mz) { out.reserve(n); for (size_t i = 0; i < n; ++i) out.push_back({(int)mz[i].pos, mz[i].val, mz[i].strand ? 1 : 0}); mz_free(mz, n); }
    return out;
}
// CSR index: a hash's hits are contiguous in hits[off[s]..off[s+1]) (sorted by sid,pos for binary search).
struct Index {
    adh::map<u64, uint32_t> id; vector<uint32_t> off; vector<Hit> hits; vector<Segment> segments;
    const Hit *range(u64 h, int &n) const {
        auto it = id.find(h);
        if (it == id.end()) { n = 0; return nullptr; }
        n = (int)(off[it->second + 1] - off[it->second]); return hits.data() + off[it->second];
    }
};
static Index build_index(const vector<std::pair<string, string>> &refs, int k, int w) {
    Index idx; vector<std::pair<u64, Hit>> raw;                    // (hash, hit) in reference order
    for (size_t sid = 0; sid < refs.size(); ++sid) {
        idx.segments.push_back({refs[sid].first, (int)refs[sid].second.size()});
        for (const auto &e : sketch(refs[sid].second, k, w)) raw.push_back({e.h, make_hit((int)sid, e.pos, e.strand)});
    }
    idx.id.reserve(raw.size());
    vector<uint32_t> slot(raw.size()), cnt;                        // dense slot per entry; hits per slot
    for (size_t i = 0; i < raw.size(); ++i) {
        uint32_t s = idx.id.try_emplace(raw[i].first, (uint32_t)cnt.size()).first->second;
        if (s == (uint32_t)cnt.size()) cnt.push_back(0);
        ++cnt[s]; slot[i] = s;
    }
    idx.off.assign(cnt.size() + 1, 0);
    for (size_t s = 0; s < cnt.size(); ++s) idx.off[s + 1] = idx.off[s] + cnt[s];  // prefix sums -> offsets
    idx.hits.resize(raw.size());
    vector<uint32_t> cur(idx.off.begin(), idx.off.end() - 1);
    for (size_t i = 0; i < raw.size(); ++i) idx.hits[cur[slot[i]]++] = raw[i].second;  // stable placement
    return idx;
}
struct Mapping { int sid, t_start, t_end; double score; int codir; int mapq = 0; bool ok = false; };
struct Seed { const Hit *hits; int n_hits; int rstrand, mult; };  // mult = read multiplicity (weighted containment)
struct BlockScore { int matches, codir, r_min, r_max; bool pruned; };
static inline double seed_heuristic(int used, int matches, int m) { return 1.0 - double(used - matches) / m; }  // upper bound in [0,1] vs theta
static inline bool disjoint(const Mapping &a, const Mapping &b) { return a.sid != b.sid || a.t_end <= b.t_start || b.t_end <= a.t_start; }  // disjoint reference intervals (Def. 5')
static vector<Seed> gather_seeds(const vector<SkEntry> &sk, const Index &idx) {  // one seed per distinct minimizer (+read multiplicity), rarest-first (stable)
    adh::map<u64, int> at; vector<Seed> seeds; at.reserve(sk.size()); seeds.reserve(sk.size());
    for (const auto &e : sk) {
        auto r = at.try_emplace(e.h, (int)seeds.size());
        if (r.second) { int nh; const Hit *hp = idx.range(e.h, nh); seeds.push_back({hp, nh, e.strand, 1}); } else ++seeds[r.first->second].mult;
    }
    std::stable_sort(seeds.begin(), seeds.end(), [](const Seed &a, const Seed &b) { return a.n_hits < b.n_hits; });
    return seeds;
}
static vector<long long> candidate_blocks(const vector<Seed> &seeds, int m, double theta, int B) {  // rarest hits -> blocks b,b-1; keys votes desc, key asc
    int S = (int)((1.0 - theta) * m) + 1;
    adh::map<long long, int> cand;                                // key (sid<<32)|block -> votes
    for (int i = 0; i < S && i < (int)seeds.size(); ++i)
        for (int j = 0; j < seeds[i].n_hits; ++j) {
            const Hit &hit = seeds[i].hits[j]; int b = hit.pos() / B;
            ++cand[((long long)hit.sid << 32) | (unsigned)b];
            if (b > 0) ++cand[((long long)hit.sid << 32) | (unsigned)(b - 1)];
        }
    vector<std::pair<int, long long>> order; order.reserve(cand.size());
    for (const auto &kv : cand) order.push_back({kv.second, kv.first});
    std::sort(order.begin(), order.end(), [](const auto &a, const auto &b) { return a.first != b.first ? a.first > b.first : a.second < b.second; });
    vector<long long> keys; keys.reserve(order.size());
    for (const auto &o : order) keys.push_back(o.second);
    return keys;
}
// Score block [lo,hi): each seed's smallest-pos hit (binary search), prune once seed_heuristic < target.
static BlockScore score_block(const vector<Seed> &seeds, int sid, int lo, int hi, int m, double target) {
    int used = 0, matches = 0, codir = 0, r_min = -1, r_max = -1;
    for (const Seed &s : seeds) {
        used += s.mult; int a = 0, z = s.n_hits;
        while (a < z) { int mid = (a + z) >> 1; const Hit &h = s.hits[mid]; if ((int)h.sid < sid || ((int)h.sid == sid && h.pos() < lo)) a = mid + 1; else z = mid; }
        if (a < s.n_hits && (int)s.hits[a].sid == sid && s.hits[a].pos() < hi) {
            const Hit &hit = s.hits[a]; int hp = hit.pos(); matches += s.mult; codir += (hit.strand() == s.rstrand) ? 1 : -1;
            r_min = r_min < 0 ? hp : std::min(r_min, hp); r_max = std::max(r_max, hp);
        }
        if (seed_heuristic(used, matches, m) < target) return {matches, codir, r_min, r_max, true};
    }
    return {matches, codir, r_min, r_max, false};
}
static Mapping map_read(const string &seq, const Index &idx, int k, int w, double theta, double delta) {  // best block + phi-free mapq, or ok=false
    Mapping best; auto sk = sketch(seq, k, w); int m = (int)sk.size();
    if (m == 0) return best;
    int B = std::max((int)seq.size(), 1);                         // candidate blocks are read-length-wide
    vector<Seed> seeds = gather_seeds(sk, idx);
    vector<long long> blocks = candidate_blocks(seeds, m, theta, B);
    long long best_key = 0; vector<Mapping> cands;                // blocks that may be best OR a near-best alternative
    for (long long key : blocks) {
        int sid = (int)(key >> 32), b = (int)(unsigned)(key & 0xffffffffLL);
        double target = best.ok ? std::max(theta, best.score - delta) : theta;   // keep within delta of best so a disjoint 2nd-best survives
        BlockScore bs = score_block(seeds, sid, b * B, (b + 2) * B, m, target);
        if (bs.pruned) continue;
        double score = double(bs.matches) / m;
        if (score < theta) continue;
        Mapping mp{sid, bs.r_min, bs.r_max + k, score, bs.codir, 0, true};
        cands.push_back(mp);
        if (!best.ok || score > best.score || (score == best.score && key < best_key)) { best = mp; best_key = key; }
    }
    if (!best.ok) return best;
    double second = 0.0;                                          // strongest alternative DISJOINT from best (Def. 6')
    for (const Mapping &mp : cands) if (disjoint(mp, best)) second = std::max(second, mp.score);
    best.mapq = (second < best.score - delta) ? 60 : 0;
    return best;
}
static string first_token(const string &s) { size_t p = s.find_first_of(" \t"); return p == string::npos ? s : s.substr(0, p); }
static vector<std::pair<string, string>> read_fasta(const string &path) {   // (name, UPPER seq); name = first header token
    vector<std::pair<string, string>> recs; std::ifstream f(path);
    if (!f) { std::cerr << "Cannot open " << path << "\n"; std::exit(1); }
    string line, name, seq; bool have = false;
    while (std::getline(f, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        if (line[0] == '>') { if (have) recs.push_back({name, seq}); name = first_token(line.substr(1)); seq.clear(); have = true; }
        else { for (char &c : line) c = (char)std::toupper((unsigned char)c); seq += line; }
    }
    if (have) recs.push_back({name, seq});
    return recs;
}
// reads[lo,hi) -> out[lo,hi); threads write disjoint slices, no locks -> byte-identical.
static void map_chunk(const vector<std::pair<string, string>> &reads, size_t lo, size_t hi, const Index &idx, int k, int w, double theta, double delta, vector<string> &out) {
    for (size_t i = lo; i < hi; ++i) {
        Mapping best = map_read(reads[i].second, idx, k, w, theta, delta);
        if (!best.ok) continue;
        const Segment &seg = idx.segments[best.sid];
        int qlen = (int)reads[i].second.size(), nmatch = (int)std::rint(best.score * qlen);    // round-half-even like Python
        std::ostringstream os;
        os << reads[i].first << '\t' << qlen << '\t' << 0 << '\t' << qlen << '\t' << (best.codir >= 0 ? '+' : '-') << '\t'
           << seg.name << '\t' << seg.length << '\t' << best.t_start << '\t' << best.t_end << '\t' << nmatch << '\t' << (best.t_end - best.t_start) << '\t' << best.mapq << '\n';
        out[i] = os.str();
    }
}
template <class F> static void parallel_for(size_t n, int threads, F work) {  // run work over [0,n) in `threads` stripes; <=1 inline
    if (threads <= 1 || n < 2) { work(0, n); return; }
    threads = (int)std::min<size_t>(threads, n);
    vector<std::thread> pool; size_t chunk = (n + threads - 1) / threads;
    for (int t = 0; t < threads; ++t) { size_t lo = std::min(n, (size_t)t * chunk), hi = std::min(n, lo + chunk); if (lo < hi) pool.emplace_back(work, lo, hi); }
    for (auto &th : pool) th.join();
}
int main(int argc, char **argv) {
    int k = 15, w = 11, threads = 1; double theta = 0.9, delta = 0.15; vector<string> pos;
    for (int i = 1; i < argc; ++i) {
        string a = argv[i]; auto next = [&]() { return string(argv[++i]); };
        if (a == "-k") k = std::stoi(next());
        else if (a == "-w" || a == "--window") w = std::stoi(next());
        else if (a == "-t" || a == "--theta") theta = std::stod(next());
        else if (a == "-d" || a == "--delta") delta = std::stod(next());
        else if (a == "-j" || a == "--threads") threads = std::stoi(next());
        else if (!a.empty() && a[0] != '-') pos.push_back(a);
        else { std::cerr << "Unknown option: " << a << "\n"; return 1; }
    }
    if (pos.size() < 2) { std::cerr << "Usage: minshmap ref.fa reads.fa [-k][-w][-t][-d][-j]\n"; return 1; }
    if (w % 2 == 0) { std::cerr << "window (-w) must be odd for canonical minimizers\n"; return 1; }
    Bench bench;
    Index idx = build_index(read_fasta(pos[0]), k, w);           // INDEX phase: read reference + build CSR index
    bench.mark_index();                                          // index time + resident memory held by the index
    auto reads = read_fasta(pos[1]); vector<string> out(reads.size());
    bench.start_map();
    parallel_for(reads.size(), threads, [&](size_t lo, size_t hi) { map_chunk(reads, lo, hi, idx, k, w, theta, delta, out); });  // MAP phase
    bench.mark_map();
    size_t mapped = 0; for (const string &line : out) if (!line.empty()) ++mapped;
    for (const string &line : out) std::cout << line;
    bench.report(reads.size(), mapped);                          // stderr-only; stdout PAF stays byte-identical
    return 0;
}
