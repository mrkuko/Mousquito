#pragma once

#include <chrono>
#include <cmath>
#include <cstdint>
#include <random>

namespace mq {

// An inclusive range [lo, hi], e.g. "2-4" from the command line.
struct Range {
    double lo = 0.0;
    double hi = 0.0;
};

// Seeded random numbers that come out the same on every OS.
//
// std::mt19937_64 is fully specified by the C++ standard, so its raw output is identical on
// libstdc++, libc++ and MSVC. The std:: distributions (normal_distribution etc.) are not:
// each library uses its own algorithm. That is why all distributions here are written by
// hand on top of the raw output. Values that go through std::log/std::cos may still differ
// in the last bit between C libraries, which never changes a rounded pixel in practice.
//
// Every call consumes a fixed number of raw values, whatever its arguments (chance(0) draws
// just like chance(50)). Changing one parameter therefore does not reshuffle all the other
// random decisions of a move with the same --seed.
class Rng {
public:
    explicit Rng(uint64_t seed) : engine_(seed) {}

    // Raw 64-bit value. Consumes 1.
    uint64_t nextU64() { return engine_(); }

    // Uniform in [0, 1) with 53 bits of precision. Consumes 1.
    double uniform01() { return static_cast<double>(nextU64() >> 11) * 0x1.0p-53; }

    // Uniform in [a, b). Consumes 1.
    double uniform(double a, double b) { return a + (b - a) * uniform01(); }
    double uniform(Range r) { return uniform(r.lo, r.hi); }

    // Uniform integer in [lo, hi], both ends included. Consumes 1 (rarely more: rejection
    // sampling keeps it unbiased).
    int uniformInt(int lo, int hi) {
        if (hi <= lo) {
            nextU64();
            return lo;
        }
        const uint64_t span = static_cast<uint64_t>(static_cast<int64_t>(hi) - lo) + 1;
        const uint64_t threshold = (0 - span) % span; // 2^64 mod span
        uint64_t x;
        do {
            x = nextU64();
        } while (x < threshold);
        return static_cast<int>(static_cast<int64_t>(lo) + static_cast<int64_t>(x % span));
    }

    // True with probability pct/100. Consumes 1.
    bool chance(double pct) { return uniform(0.0, 100.0) < pct; }

    // Gaussian (Box–Muller). Consumes 2.
    double normal(double mean, double sigma) {
        const double u1 = 1.0 - uniform01(); // (0, 1], so log() is finite
        const double u2 = uniform01();
        const double z = std::sqrt(-2.0 * std::log(u1)) * std::cos(kTwoPi * u2);
        return mean + sigma * z;
    }

    // Log-normal factor with mean 1, the same jitter as Pasteguise: right-skewed (occasional
    // long values, never negative) without changing the average. sigma = jitter% / 100.
    // sigma 0 returns exactly 1. Consumes 2.
    double logNormalMean1(double sigma) {
        const double n = normal(0.0, sigma);
        return sigma > 0.0 ? std::exp(n - sigma * sigma / 2.0) : 1.0;
    }

private:
    static constexpr double kTwoPi = 6.283185307179586476925286766559;

    std::mt19937_64 engine_;
};

// A seed for runs without --seed. Printed by the CLI so a run can be replayed.
inline uint64_t randomSeed() {
    std::random_device rd;
    const uint64_t hw = (static_cast<uint64_t>(rd()) << 32) ^ rd();
    const auto now = static_cast<uint64_t>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    return hw ^ (now * 0x9E3779B97F4A7C15ull);
}

} // namespace mq
