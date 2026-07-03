// Optional benchmark instrumentation for minSHmap - kept OUT of the pedagogical mapper.
//
// minshmap.cpp stays minimal: it makes a `Bench` in main(), calls mark_index() / start_map() /
// mark_map() at the phase boundaries, then report(). All the chrono timing, /proc RSS parsing
// and formatting live here so they don't clutter the teaching file.
//
// Everything is a near-no-op unless the environment variable MINSHMAP_BENCH is set; when it is,
// exactly ONE line goes to STDERR (stdout PAF is never touched, so py/cpp stay byte-identical):
//
//   [bench] index_s=.. map_s=.. reads=.. mapped=.. index_rss_mb=.. peak_rss_mb=..
//
// The real-world benchmark (realworld/11_bench_3way.py) sets MINSHMAP_BENCH and regex-parses it.
#ifndef MINSHMAP_BENCH_HPP
#define MINSHMAP_BENCH_HPP
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

struct Bench {
    using clk = std::chrono::steady_clock;
    bool on = std::getenv("MINSHMAP_BENCH") != nullptr;   // off -> mark/report do nothing
    clk::time_point lap = clk::now();                     // start of the current phase
    double index_s = 0, map_s = 0, index_rss_mb = 0;

    // current VmRSS / peak VmHWM in MB from /proc/self/status (Linux); 0 elsewhere.
    static void rss_mb(double &cur_mb, double &peak_mb) {
        cur_mb = peak_mb = 0.0; std::ifstream f("/proc/self/status"); std::string line;
        while (std::getline(f, line)) {
            if (line.rfind("VmRSS:", 0) == 0) cur_mb = std::stod(line.substr(6)) / 1024.0;
            else if (line.rfind("VmHWM:", 0) == 0) peak_mb = std::stod(line.substr(6)) / 1024.0;
        }
    }
    void mark_index() {                                   // close INDEX phase: duration + resident memory of the index
        if (!on) return;
        auto now = clk::now(); index_s = std::chrono::duration<double>(now - lap).count();
        double peak; rss_mb(index_rss_mb, peak); lap = now;
    }
    void start_map() { if (on) lap = clk::now(); }        // begin MAP phase (excludes reading the reads file)
    void mark_map()  { if (on) map_s = std::chrono::duration<double>(clk::now() - lap).count(); }
    void report(size_t reads, size_t mapped) {
        if (!on) return;
        double cur = 0, peak = 0; rss_mb(cur, peak);
        std::cerr << "[bench] index_s=" << index_s << " map_s=" << map_s
                  << " reads=" << reads << " mapped=" << mapped
                  << " index_rss_mb=" << index_rss_mb << " peak_rss_mb=" << peak << "\n";
    }
};
#endif
