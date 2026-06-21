// minSHmap (C++ port) - a minimalistic sketch-based read mapper.
//
// A line-for-line equivalent of minshmap.py, kept in the same spirit as `minSH`
// and `shmap`: simple to understand, with the *seed heuristic* (the "SH") at its
// core. It exists so the Python and C++ implementations can be compared on the
// same data and the same algorithm.
//
//   1. SKETCH - FracMinHash: keep the small fraction `hfrac` of k-mer hashes.
//   2. INDEX  - reference hash -> [(segment, position, strand), ...].
//   3. MAP    - for each read, find the reference window sharing the most
//               k-mers, pruning hopeless windows with
//                   sh = 1 - (seeds_used - matches) / m
//               (an upper bound on the containment a window can still reach).
//
// Three interchangeable sketchers (select with --hash): naive | poly | nthash.
//
// Build:  g++ -O3 -std=c++17 -march=native -o minshmap minshmap.cpp
// Usage:  ./minshmap ref.fa reads.fa [-k 15] [-r 0.05] [-t 0.9] [--hash nthash]
//         ./minshmap --demo [--hash poly]
//         ./minshmap ref.fa reads.fa --report   # tsv timing+accuracy (truth in headers)

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using u64 = uint64_t;
using std::string;
using std::vector;

static const u64 MASK64 = ~u64(0);

// --------------------------------------------------------------------------- //
// Small helpers
// --------------------------------------------------------------------------- //

static inline char comp(char c) {
    switch (c) { case 'A': return 'T'; case 'C': return 'G';
                 case 'G': return 'C'; case 'T': return 'A'; default: return 'N'; }
}
static string revcomp(const string &s) {
    string r(s.size(), 'N');
    for (size_t i = 0; i < s.size(); ++i) r[s.size() - 1 - i] = comp(s[i]);
    return r;
}
static inline int base_val(char c) {
    switch (c) { case 'A': return 0; case 'C': return 1;
                 case 'G': return 2; case 'T': return 3; default: return 0; }
}
static inline u64 thr64(double hfrac) {
    return hfrac >= 1.0 ? MASK64 : (u64)(hfrac * 18446744073709551616.0);  // hfrac * 2^64
}

// One kept k-mer of a sketch.
struct SkEntry { int pos; u64 h; int strand; };  // strand: 0 forward, 1 rev-comp

// --------------------------------------------------------------------------- //
// Sketchers (FracMinHash)
// --------------------------------------------------------------------------- //

// ---- 1. Naive: hash each k-mer string from scratch.                 O(len*k)
static u64 fnv1a(const char *s, int n) {
    u64 h = 1469598103934665603ULL;
    for (int i = 0; i < n; ++i) { h ^= (unsigned char)s[i]; h *= 1099511628211ULL; }
    return h;
}
static vector<SkEntry> sketch_naive(const string &seq, int k, double hfrac) {
    vector<SkEntry> out;
    u64 thr = thr64(hfrac);
    int n = (int)seq.size();
    string rc(k, 'N');
    for (int i = 0; i + k <= n; ++i) {
        for (int j = 0; j < k; ++j) rc[k - 1 - j] = comp(seq[i + j]);
        u64 h_fw = fnv1a(seq.data() + i, k);
        u64 h_rc = fnv1a(rc.data(), k);
        u64 h = std::min(h_fw, h_rc);  // canonical
        if (h <= thr) out.push_back({i, h, h_fw <= h_rc ? 0 : 1});
    }
    return out;
}

// ---- 2. Polynomial (Rabin-Karp) rolling hash.                          O(len)
static const u64 POLY_M = (u64(1) << 61) - 1;   // Mersenne prime
static const u64 POLY_B = 0x9E3779B1ULL;
static inline u64 mulmod(u64 a, u64 b) { return (u64)((__uint128_t)a * b % POLY_M); }
static u64 powmod(u64 a, u64 e) { u64 r = 1; a %= POLY_M; while (e) { if (e & 1) r = mulmod(r, a); a = mulmod(a, a); e >>= 1; } return r; }
static vector<SkEntry> sketch_poly(const string &seq, int k, double hfrac) {
    vector<SkEntry> out;
    int n = (int)seq.size();
    if (n < k) return out;
    u64 thr = (u64)(hfrac * (double)POLY_M);
    u64 invB = powmod(POLY_B, POLY_M - 2);
    u64 bk1 = powmod(POLY_B, k - 1);

    u64 h_fw = 0, h_rc = 0, Bt = 1;
    for (int t = 0; t < k; ++t) {
        int v = base_val(seq[t]);
        h_fw = (mulmod(h_fw, POLY_B) + v) % POLY_M;
        h_rc = (h_rc + mulmod((u64)(3 - v), Bt)) % POLY_M;
        Bt = mulmod(Bt, POLY_B);
    }
    for (int i = 0; i + k <= n; ++i) {
        u64 h = std::min(h_fw, h_rc);
        if (h <= thr) out.push_back({i, h, h_fw <= h_rc ? 0 : 1});
        if (i + k < n) {
            int out_v = base_val(seq[i]), in_v = base_val(seq[i + k]);
            h_fw = (mulmod((h_fw + POLY_M - mulmod((u64)out_v, bk1)) % POLY_M, POLY_B) + in_v) % POLY_M;
            u64 rc_drop = (h_rc + POLY_M - (u64)(3 - out_v)) % POLY_M;
            h_rc = (mulmod(rc_drop, invB) + mulmod((u64)(3 - in_v), bk1)) % POLY_M;
        }
    }
    return out;
}

// ---- 3. ntHash-style rotation rolling hash (as in shmap).              O(len)
static u64 LUT_FW[256], LUT_RC[256];
static void init_lut() {
    LUT_FW['A'] = 0x3C8BFBB395C60474ULL; LUT_FW['C'] = 0x3193C18562A02B4CULL;
    LUT_FW['G'] = 0x20323ED082572324ULL; LUT_FW['T'] = 0x295549F54BE24456ULL;
    LUT_RC['A'] = LUT_FW['T']; LUT_RC['C'] = LUT_FW['G'];
    LUT_RC['G'] = LUT_FW['C']; LUT_RC['T'] = LUT_FW['A'];
}
static inline u64 rotl(u64 x, int r) { r &= 63; return r ? (x << r) | (x >> (64 - r)) : x; }
static inline u64 rotr(u64 x, int r) { r &= 63; return r ? (x >> r) | (x << (64 - r)) : x; }
static vector<SkEntry> sketch_nthash(const string &seq, int k, double hfrac) {
    vector<SkEntry> out;
    int n = (int)seq.size();
    if (n < k) return out;
    u64 thr = thr64(hfrac);

    u64 h_fw = 0, h_rc = 0;
    for (int t = 0; t < k; ++t) {
        h_fw ^= rotl(LUT_FW[(unsigned char)seq[t]], k - 1 - t);
        h_rc ^= rotl(LUT_RC[(unsigned char)seq[t]], t);
    }
    for (int i = 0; i + k <= n; ++i) {
        u64 h = std::min(h_fw, h_rc);  // canonical: min of the two strands
        if (h <= thr) out.push_back({i, h, h_fw > h_rc ? 1 : 0});
        if (i + k < n) {
            unsigned char o = seq[i], in = seq[i + k];
            h_fw = rotl(h_fw, 1) ^ rotl(LUT_FW[o], k) ^ LUT_FW[in];
            h_rc = rotr(h_rc, 1) ^ rotr(LUT_RC[o], 1) ^ rotl(LUT_RC[in], k - 1);
        }
    }
    return out;
}

using Sketcher = vector<SkEntry> (*)(const string &, int, double);
static Sketcher pick_sketcher(const string &name) {
    if (name == "naive") return sketch_naive;
    if (name == "poly") return sketch_poly;
    return sketch_nthash;
}

// --------------------------------------------------------------------------- //
// Index + mapping
// --------------------------------------------------------------------------- //

struct Hit { int sid; int pos; int strand; };
struct Segment { string name; int length; };
struct Index { std::unordered_map<u64, vector<Hit>> h2hits; vector<Segment> segments; };

static Index build_index(const vector<std::pair<string, string>> &refs, int k, double hfrac, Sketcher sketch) {
    Index idx;
    for (size_t sid = 0; sid < refs.size(); ++sid) {
        idx.segments.push_back({refs[sid].first, (int)refs[sid].second.size()});
        for (const auto &e : sketch(refs[sid].second, k, hfrac))
            idx.h2hits[e.h].push_back({(int)sid, e.pos, e.strand});
    }
    return idx;
}

struct Mapping { int sid; int t_start; int t_end; double score; int codir; bool ok = false; };

// Containment of the read in window [lo,hi) of segment `sid`, with seed-heuristic
// pruning. Returns ok=false as soon as the heuristic falls below `thr`.
struct ShResult { bool ok; double score; int codir; int rmin; int rmax; };
static ShResult containment_with_sh(
        const vector<std::tuple<int, int, u64, const vector<Hit> *, int>> &seeds,
        int m, int sid, int lo, int hi, double thr) {
    int seeds_used = 0, matches = 0, codir = 0, rmin = INT32_MAX, rmax = -1;
    for (const auto &s : seeds) {
        int occ_in_p = std::get<1>(s);
        int rstrand = std::get<4>(s);
        const vector<Hit> *hits = std::get<3>(s);
        seeds_used += occ_in_p;
        int matched = 0;
        if (hits) {
            for (const Hit &hit : *hits) {
                if (hit.sid == sid && lo <= hit.pos && hit.pos < hi) {
                    ++matched;
                    codir += (hit.strand == rstrand) ? 1 : -1;
                    rmin = std::min(rmin, hit.pos);
                    rmax = std::max(rmax, hit.pos);
                }
            }
        }
        matches += std::min(matched, occ_in_p);
        if (1.0 - double(seeds_used - matches) / m < thr) return {false, 0, 0, 0, 0};
    }
    return {true, double(matches) / m, codir, rmin, rmax};
}

static std::pair<Mapping, Mapping> map_read(const string &seq, const Index &idx, int k,
                                            double hfrac, double theta, double min_diff,
                                            Sketcher sketch, int max_matches) {
    Mapping none, none2;
    auto sk = sketch(seq, k, hfrac);
    int m = (int)sk.size();
    if (m == 0) return {none, none2};
    int halflen = std::max((int)seq.size(), 1);

    // Group read k-mers by hash: multiplicity and (first-seen) strand.
    std::unordered_map<u64, int> occ, rstrand;
    for (const auto &e : sk) { occ[e.h]++; rstrand.emplace(e.h, e.strand); }

    // One seed per distinct read k-mer, annotated with rarity in the index.
    // tuple = (n_hits, occ_in_p, hash, hits*, read_strand)
    vector<std::tuple<int, int, u64, const vector<Hit> *, int>> seeds;
    seeds.reserve(occ.size());
    for (const auto &kv : occ) {
        auto it = idx.h2hits.find(kv.first);
        const vector<Hit> *hits = (it == idx.h2hits.end()) ? nullptr : &it->second;
        int n_hits = hits ? (int)hits->size() : 0;
        // Drop over-frequent (uninformative) k-mers; recompute informative m.
        if (max_matches > 0 && n_hits > max_matches) continue;
        seeds.emplace_back(n_hits, kv.second, kv.first, hits, rstrand[kv.first]);
    }
    if (max_matches > 0) {
        m = 0;
        for (const auto &s : seeds) m += std::get<1>(s);
        if (m == 0) return {none, none2};
    }
    std::sort(seeds.begin(), seeds.end(),
              [](const auto &a, const auto &b) { return std::get<0>(a) < std::get<0>(b); });

    // Seed candidate windows using the rarest seeds (overlapping b and b-1).
    double theta2 = theta - min_diff;
    int S = (int)((1.0 - theta2) * m) + 1;
    std::unordered_set<long long> candidates;
    int used = 0;
    for (const auto &s : seeds) {
        if (used >= S) break;
        int n_hits = std::get<0>(s);
        const vector<Hit> *hits = std::get<3>(s);
        if (n_hits == 0) continue;
        for (const Hit &hit : *hits) {
            int b = hit.pos / halflen;
            candidates.insert(((long long)hit.sid << 32) | (unsigned)b);
            if (b > 0) candidates.insert(((long long)hit.sid << 32) | (unsigned)(b - 1));
        }
        ++used;
    }

    // Refine each candidate window, pruning with the seed heuristic.
    vector<Mapping> found;
    for (long long key : candidates) {
        int sid = (int)(key >> 32);
        int b = (int)(unsigned)(key & 0xffffffff);
        int lo = b * halflen, hi = (b + 2) * halflen;
        ShResult r = containment_with_sh(seeds, m, sid, lo, hi, theta);
        if (!r.ok || r.score < theta) continue;
        found.push_back({sid, r.rmin, r.rmax + k, r.score, r.codir, true});
    }
    if (found.empty()) return {none, none2};
    std::sort(found.begin(), found.end(),
              [](const Mapping &a, const Mapping &b) { return a.score > b.score; });
    Mapping best = found[0], second;
    for (size_t i = 1; i < found.size(); ++i) {
        const Mapping &mp = found[i];
        if (mp.sid != best.sid || mp.t_start >= best.t_end || mp.t_end <= best.t_start) { second = mp; break; }
    }
    return {best, second};
}

static int mapq(const Mapping &best, const Mapping &second, double theta, double min_diff) {
    double score2 = second.ok ? second.score : (theta - min_diff);
    double frac = 1.0 - score2 / best.score;
    return frac > min_diff ? 60 : 0;
}

// --------------------------------------------------------------------------- //
// I/O
// --------------------------------------------------------------------------- //

// Reads FASTA keeping the full header line (truth lives after the read name).
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
            name = line.substr(1); seq.clear(); have = true;
        } else {
            for (char &c : line) c = (char)std::toupper((unsigned char)c);
            seq += line;
        }
    }
    if (have) recs.push_back({name, seq});
    return recs;
}

static string first_token(const string &s) {
    size_t p = s.find_first_of(" \t");
    return p == string::npos ? s : s.substr(0, p);
}

static void print_paf(const string &qname, int qlen, const Mapping &best, const Mapping &second,
                      const Index &idx, double theta, double min_diff) {
    const Segment &seg = idx.segments[best.sid];
    char strand = best.codir >= 0 ? '+' : '-';
    int nmatch = (int)(best.score * qlen + 0.5);
    std::cout << qname << '\t' << qlen << '\t' << 0 << '\t' << qlen << '\t' << strand << '\t'
              << seg.name << '\t' << seg.length << '\t' << best.t_start << '\t' << best.t_end << '\t'
              << nmatch << '\t' << (best.t_end - best.t_start) << '\t'
              << mapq(best, second, theta, min_diff) << '\n';
}

// --------------------------------------------------------------------------- //
// Demo / report / main
// --------------------------------------------------------------------------- //

static string random_seq(int n, std::mt19937 &rng) {
    static const char *N = "ACGT";
    string s(n, 'A');
    for (int i = 0; i < n; ++i) s[i] = N[rng() & 3];
    return s;
}
static string mutate(const string &s, double e, std::mt19937 &rng) {
    static const char *N = "ACGT";
    std::uniform_real_distribution<double> U(0, 1);
    string out;
    for (char c : s) {
        if (U(rng) > e) out += c;
        else if (U(rng) < 1.0 / 3) out += N[rng() & 3];
        else if (U(rng) < 1.0 / 2) { out += N[rng() & 3]; out += c; }
        // else: deletion
    }
    return out;
}

static int run_demo(Sketcher sketch) {
    // Demo-specific parameters (match minshmap.py demo): a read k-mer survives a
    // mutated read only if all k bases are error-free, so expected containment
    // ~= (1-e)^k must exceed theta. Here e=4%: (1-0.04)^11 ~ 0.64 > theta=0.5.
    int k = 11; double hfrac = 0.2, theta = 0.5, min_diff = 0.02, error = 0.04;
    std::mt19937 rng(1);
    string ref = random_seq(5000, rng);
    Index idx = build_index({{"chr_demo", ref}}, k, hfrac, sketch);
    std::uniform_int_distribution<int> startD(0, (int)ref.size() - 400);
    std::uniform_real_distribution<double> U(0, 1);
    int correct = 0, total = 20;
    for (int i = 0; i < total; ++i) {
        int start = startD(rng);
        string read = mutate(ref.substr(start, 300), error, rng);
        if (U(rng) < 0.5) read = revcomp(read);
        auto [best, second] = map_read(read, idx, k, hfrac, theta, min_diff, sketch, 0);
        bool ok = best.ok && best.t_start <= start + 150 && start + 150 <= best.t_end;
        correct += ok;
        std::printf("read%02d true~%5d  ->  %s  %s\n", i, start,
                    best.ok ? (std::to_string(best.t_start) + "-" + std::to_string(best.t_end) +
                               " score=" + std::to_string(best.score)).c_str()
                            : "UNMAPPED",
                    ok ? "OK" : "MISS");
    }
    std::printf("\nCorrectly placed %d/%d reads.\n", correct, total);
    return 0;
}

// Parse "segm=..", "pos=..", "strand=.." from a read header.
static void parse_truth(const string &header, string &segm, int &pos) {
    segm.clear(); pos = -1;
    std::istringstream ss(header);
    string tok;
    while (ss >> tok) {
        auto eq = tok.find('=');
        if (eq == string::npos) continue;
        string key = tok.substr(0, eq), val = tok.substr(eq + 1);
        if (key == "segm") segm = val;
        else if (key == "pos") pos = std::stoi(val);
    }
}

static int run_report(const string &ref_path, const string &reads_path, int k, double hfrac,
                      double theta, double min_diff, const string &hash_name, Sketcher sketch,
                      int max_matches) {
    auto refs = read_fasta(ref_path);
    auto reads = read_fasta(reads_path);
    long ref_bp = 0;
    for (auto &r : refs) ref_bp += (long)r.second.size();
    std::unordered_map<string, int> seg_id;

    using clk = std::chrono::steady_clock;
    auto t0 = clk::now();
    Index idx = build_index(refs, k, hfrac, sketch);
    double t_index = std::chrono::duration<double>(clk::now() - t0).count();
    for (size_t i = 0; i < idx.segments.size(); ++i) seg_id[idx.segments[i].name] = (int)i;

    int mapped = 0, correct = 0;
    t0 = clk::now();
    for (auto &rd : reads) {
        auto [best, second] = map_read(rd.second, idx, k, hfrac, theta, min_diff, sketch, max_matches);
        if (!best.ok) continue;
        ++mapped;
        string segm; int pos;
        parse_truth(rd.first, segm, pos);
        auto it = seg_id.find(segm);
        int mid = (int)rd.second.size() / 2;
        if (it != seg_id.end() && best.sid == it->second &&
            best.t_start <= pos + mid && pos + mid <= best.t_end)
            ++correct;
    }
    double t_map = std::chrono::duration<double>(clk::now() - t0).count();

    int n = (int)reads.size();
    // tsv: hash index_sec ref_bp_per_s map_sec reads_per_s mapped_frac accuracy
    std::printf("%s\t%.4f\t%ld\t%.4f\t%ld\t%.3f\t%.3f\n",
                hash_name.c_str(), t_index,
                t_index > 0 ? (long)(ref_bp / t_index) : 0,
                t_map, t_map > 0 ? (long)(n / t_map) : 0,
                n ? double(mapped) / n : 0.0, n ? double(correct) / n : 0.0);
    return 0;
}

int main(int argc, char **argv) {
    init_lut();
    string ref, reads, hash_name = "nthash";
    int k = 15;
    double hfrac = 0.05, theta = 0.9, min_diff = 0.02;
    int max_matches = 0;
    bool demo = false, report = false;
    vector<string> pos;

    for (int i = 1; i < argc; ++i) {
        string a = argv[i];
        auto next = [&]() { return string(argv[++i]); };
        if (a == "--demo") demo = true;
        else if (a == "--report") report = true;
        else if (a == "--hash") hash_name = next();
        else if (a == "-k") k = std::stoi(next());
        else if (a == "-r" || a == "--hfrac") hfrac = std::stod(next());
        else if (a == "-t" || a == "--theta") theta = std::stod(next());
        else if (a == "-d" || a == "--min-diff") min_diff = std::stod(next());
        else if (a == "-M" || a == "--max-matches") max_matches = std::stoi(next());
        else if (!a.empty() && a[0] != '-') pos.push_back(a);
        else { std::cerr << "Unknown option: " << a << "\n"; return 1; }
    }
    Sketcher sketch = pick_sketcher(hash_name);

    if (demo) return run_demo(sketch);
    if (pos.size() < 2) { std::cerr << "Usage: minshmap ref.fa reads.fa [options] | --demo\n"; return 1; }
    ref = pos[0]; reads = pos[1];
    if (report) return run_report(ref, reads, k, hfrac, theta, min_diff, hash_name, sketch, max_matches);

    Index idx = build_index(read_fasta(ref), k, hfrac, sketch);
    for (auto &rd : read_fasta(reads)) {
        auto [best, second] = map_read(rd.second, idx, k, hfrac, theta, min_diff, sketch, max_matches);
        if (best.ok) print_paf(first_token(rd.first), (int)rd.second.size(), best, second, idx, theta, min_diff);
    }
    return 0;
}
