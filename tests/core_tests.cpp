#include <cmath>
#include <cstdio>
#include <cstring>

#include "geometry.hpp"
#include "motion.hpp"
#include "random.hpp"
#include "trajectory.hpp"
#include "version.hpp"

static int failures = 0;
#define CHECK(cond)                                                          \
    do {                                                                     \
        if (!(cond)) {                                                       \
            std::fprintf(stderr, "FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond); \
            ++failures;                                                      \
        }                                                                    \
    } while (0)

static bool approx(double a, double b, double eps) { return std::fabs(a - b) <= eps; }
static bool approx(mq::Vec2 a, mq::Vec2 b, double eps) { return approx(a.x, b.x, eps) && approx(a.y, b.y, eps); }

static void testVersion() {
    const char* v = mq::version();
    CHECK(v != nullptr);
    CHECK(std::strlen(v) >= 5); // "X.Y.Z"
}

static void testVectors() {
    using mq::Vec2;
    CHECK(approx(Vec2(3, 4).length(), 5.0, 1e-12));
    CHECK(approx(mq::normalized(Vec2(0, 7)), Vec2(0, 1), 1e-12));
    CHECK(approx(mq::normalized(Vec2()), Vec2(), 0.0));

    // The normal is a unit vector at 90° to the direction, and zero for a zero direction.
    const Vec2 dir(5, -2);
    const Vec2 n = mq::unitNormal(dir);
    CHECK(approx(n.length(), 1.0, 1e-12));
    CHECK(approx(mq::dot(n, dir), 0.0, 1e-12));
    CHECK(approx(mq::unitNormal(Vec2()), Vec2(), 0.0));

    CHECK(mq::roundToPoint(Vec2(2.5, -2.5)) == (mq::Point{3, -3}));
    CHECK(mq::roundToPoint(Vec2(-0.4, 0.4)) == (mq::Point{0, 0}));
}

static void testRect() {
    // A desktop with a monitor left of the primary one.
    const mq::Rect r{-1920, 0, 3840, 1080};
    CHECK(r.contains({-1920, 0}));
    CHECK(r.contains({1919, 1079}));
    CHECK(!r.contains({1920, 0}));
    CHECK(!r.contains({0, -1}));
    CHECK(r.clamp({-5000, 500}) == (mq::Point{-1920, 500}));
    CHECK(r.clamp({5000, 5000}) == (mq::Point{1919, 1079}));
    CHECK(r.clamp({10, 10}) == (mq::Point{10, 10}));

    const mq::Rect empty{};
    CHECK(empty.empty());
    CHECK(!empty.contains({0, 0}));
    CHECK(empty.clamp({7, 8}) == (mq::Point{7, 8}));
}

static void testBezier() {
    using mq::Vec2;
    const mq::CubicBezier c{{10, 20}, {50, -40}, {200, 300}, {400, 100}};
    CHECK(approx(c.at(0.0), c.p0, 1e-12));
    CHECK(approx(c.at(1.0), c.p3, 1e-12));

    // The derivative at the ends points along the control polygon: 3·(P1−P0) and 3·(P3−P2).
    CHECK(approx(c.derivative(0.0), (c.p1 - c.p0) * 3.0, 1e-9));
    CHECK(approx(c.derivative(1.0), (c.p3 - c.p2) * 3.0, 1e-9));

    // The derivative matches a finite difference in the middle.
    const double h = 1e-6;
    const Vec2 fd = (c.at(0.4 + h) - c.at(0.4 - h)) / (2 * h);
    CHECK(approx(c.derivative(0.4), fd, 1e-4));
}

static void testArcLength() {
    using mq::Vec2;

    // Straight line with unevenly placed control points: length is the distance between
    // the ends, even though equal t steps are unequal in space here.
    const mq::CubicBezier line{{0, 0}, {30, 40}, {270, 360}, {300, 400}};
    const mq::ArcLengthTable lt(line);
    CHECK(approx(lt.length(), 500.0, 1e-9));

    // Quarter circle of radius 100 (standard Bézier approximation): length ≈ π/2 · 100.
    const double k = 0.5522847498 * 100;
    const mq::CubicBezier arc{{100, 0}, {100, k}, {k, 100}, {0, 100}};
    const mq::ArcLengthTable at(arc);
    CHECK(approx(at.length(), 3.14159265358979 / 2 * 100, 0.05));

    // Clamping at both ends.
    CHECK(at.tAtLength(-5) == 0.0);
    CHECK(at.tAtLength(0) == 0.0);
    CHECK(at.tAtLength(at.length()) == 1.0);
    CHECK(at.tAtLength(at.length() + 5) == 1.0);
    CHECK(at.lengthAtT(-1) == 0.0);
    CHECK(approx(at.lengthAtT(2), at.length(), 1e-12));

    // The inversion is monotonic and round-trips t → s → t.
    double prevT = -1.0;
    bool monotonic = true, roundTrip = true;
    for (int i = 0; i <= 1000; ++i) {
        const double s = at.length() * i / 1000.0;
        const double t = at.tAtLength(s);
        if (t < prevT) monotonic = false;
        prevT = t;
        if (!approx(at.lengthAtT(t), s, 1e-9)) roundTrip = false;
    }
    CHECK(monotonic);
    CHECK(roundTrip);
    for (int i = 0; i <= 100; ++i) {
        const double t = i / 100.0;
        if (!approx(at.tAtLength(at.lengthAtT(t)), t, 1e-9)) roundTrip = false;
    }
    CHECK(roundTrip);

    // The point of the table: equal steps in s give (nearly) equal steps in space,
    // on a curve where equal steps in t clearly do not.
    const mq::ArcLengthTable ut(line);
    const int steps = 100;
    const double expected = ut.length() / steps;
    double worstS = 0.0, worstT = 0.0;
    for (int i = 0; i < steps; ++i) {
        const double dS = mq::distance(line.at(ut.tAtLength(expected * i)), line.at(ut.tAtLength(expected * (i + 1))));
        const double dT = mq::distance(line.at(double(i) / steps), line.at(double(i + 1) / steps));
        worstS = std::fmax(worstS, std::fabs(dS - expected) / expected);
        worstT = std::fmax(worstT, std::fabs(dT - expected) / expected);
    }
    CHECK(worstS < 0.01);
    CHECK(worstT > 0.1);

    // A curve that does not move has length 0 and never divides by zero.
    const mq::CubicBezier dot{{5, 5}, {5, 5}, {5, 5}, {5, 5}};
    const mq::ArcLengthTable dt(dot);
    CHECK(dt.length() == 0.0);
    CHECK(dt.tAtLength(0) == 0.0);
    CHECK(dt.tAtLength(10) == 0.0);
    CHECK(dt.lengthAtT(0.5) == 0.0);

    // A tiny table still works.
    const mq::ArcLengthTable tiny(line, 2);
    CHECK(approx(tiny.length(), 500.0, 1e-9));
    CHECK(approx(tiny.tAtLength(250), 0.5, 1e-12));
}

static void testRandom() {
    // The engine itself is the standard one: the C++ standard fixes the 10000th output of
    // a default-seeded mt19937_64.
    {
        mq::Rng r(5489);
        uint64_t x = 0;
        for (int i = 0; i < 10000; ++i) x = r.nextU64();
        CHECK(x == 9981545732273789042ull);
    }

    // Golden values for seed 7. If these fail on one OS only, --seed is not portable there.
    {
        mq::Rng r(7);
        CHECK(r.nextU64() == 13915952638675311015ull);
        CHECK(r.uniform01() == 0.94930120289264419);
        CHECK(r.uniformInt(-100, 100) == 5);
        CHECK(approx(r.normal(0, 1), 1.3315647985499977, 1e-12));
        CHECK(r.chance(50));
    }

    // Same seed, same sequence. Different seed, different sequence.
    {
        mq::Rng a(42), b(42), c(43);
        bool same = true, differs = false;
        for (int i = 0; i < 100; ++i) {
            const double va = a.normal(0, 1), vb = b.normal(0, 1), vc = c.normal(0, 1);
            if (va != vb) same = false;
            if (va != vc) differs = true;
        }
        CHECK(same);
        CHECK(differs);
    }

    const int n = 200000;

    // uniform01 stays in [0, 1) and averages 0.5. uniform(a, b) stays in [a, b).
    {
        mq::Rng r(1);
        double sum = 0, lo = 1, hi = 0;
        bool inside = true;
        for (int i = 0; i < n; ++i) {
            const double u = r.uniform01();
            sum += u;
            lo = std::fmin(lo, u);
            hi = std::fmax(hi, u);
            const double v = r.uniform(mq::Range{-3, 5});
            if (v < -3 || v >= 5) inside = false;
        }
        CHECK(lo >= 0.0 && hi < 1.0);
        CHECK(approx(sum / n, 0.5, 0.005));
        CHECK(inside);
    }

    // uniformInt hits both ends and nothing outside, also for negative and one-value ranges.
    {
        mq::Rng r(2);
        int counts[7] = {};
        bool inside = true;
        for (int i = 0; i < 70000; ++i) {
            const int v = r.uniformInt(-3, 3);
            if (v < -3 || v > 3) inside = false;
            else ++counts[v + 3];
        }
        CHECK(inside);
        bool even = true;
        for (int c : counts)
            if (c < 9000 || c > 11000) even = false;
        CHECK(even);
        CHECK(r.uniformInt(4, 4) == 4);
        CHECK(r.uniformInt(9, 2) == 9); // reversed range: lo
    }

    // normal has the requested mean and standard deviation.
    {
        mq::Rng r(3);
        double sum = 0, sq = 0;
        for (int i = 0; i < n; ++i) {
            const double v = r.normal(10, 2);
            sum += v;
            sq += v * v;
        }
        const double mean = sum / n;
        const double sd = std::sqrt(sq / n - mean * mean);
        CHECK(approx(mean, 10.0, 0.02));
        CHECK(approx(sd, 2.0, 0.02));
    }

    // logNormalMean1 averages 1, is always positive, and is exactly 1 with sigma 0.
    {
        mq::Rng r(4);
        double sum = 0, lo = 1e9, hi = 0;
        for (int i = 0; i < n; ++i) {
            const double v = r.logNormalMean1(0.3);
            sum += v;
            lo = std::fmin(lo, v);
            hi = std::fmax(hi, v);
        }
        CHECK(approx(sum / n, 1.0, 0.005));
        CHECK(lo > 0.0);
        CHECK(hi > 2.0); // right-skewed: occasional long values
        CHECK(r.logNormalMean1(0.0) == 1.0);
    }

    // chance: 0 never, 100 always, 25 about a quarter of the time.
    {
        mq::Rng r(5);
        int never = 0, always = 0, quarter = 0;
        for (int i = 0; i < n; ++i) {
            if (r.chance(0)) ++never;
            if (r.chance(100)) ++always;
            if (r.chance(25)) ++quarter;
        }
        CHECK(never == 0);
        CHECK(always == n);
        CHECK(approx(quarter / double(n), 0.25, 0.005));
    }

    // Each call consumes the same amount whatever its arguments, so changing one parameter
    // does not shift the random decisions that come after it.
    {
        mq::Rng a(6), b(6);
        a.chance(0);
        b.chance(50);
        a.logNormalMean1(0.0);
        b.logNormalMean1(0.4);
        a.uniformInt(1, 1);
        b.uniformInt(1, 6);
        CHECK(a.nextU64() == b.nextU64());
    }

    CHECK(mq::randomSeed() != mq::randomSeed());
}

static void testTrajectory() {
    // Empty and single-sample trajectories have no duration, length or speed.
    {
        mq::Trajectory t;
        CHECK(t.empty());
        CHECK(t.durationMs() == 0.0);
        CHECK(t.length() == 0.0);
        CHECK(t.straightness() == 1.0);
        CHECK(t.peakSpeed() == 0.0);
        CHECK(t.peakSpeedTimeMs() == 0.0);
        t.append({0, 5, 5});
        CHECK(t.size() == 1);
        CHECK(t.durationMs() == 0.0);
        CHECK(t.length() == 0.0);
        CHECK(t.peakSpeed() == 0.0);
    }

    // An L-shaped path: 30 px right in 30 ms, then 40 px down in 10 ms, then a pause.
    //   (0,0) -> (10,0) -> (20,0) -> (30,0) -> (30,40) -> (30,40)
    //   t = 0     10        20        30        40        50
    {
        mq::Trajectory t;
        t.append({0, 0, 0});
        t.append({10, 10, 0});
        t.append({20, 20, 0});
        t.append({30, 30, 0});
        t.append({40, 30, 40});
        t.append({50, 30, 40});
        CHECK(t.size() == 6);
        CHECK(t.front().point() == (mq::Point{0, 0}));
        CHECK(t.back().point() == (mq::Point{30, 40}));
        CHECK(approx(t.durationMs(), 50.0, 1e-12));
        CHECK(approx(t.length(), 70.0, 1e-12));
        CHECK(approx(t.straightness(), 50.0 / 70.0, 1e-12)); // 30-40-50 triangle
        CHECK(approx(t.peakSpeed(), 4.0, 1e-12));            // 40 px in 10 ms
        CHECK(approx(t.peakSpeedTimeMs(), 35.0, 1e-12));     // middle of that step
    }

    // A diagonal step counts its true length, and times are relative to the first sample.
    {
        mq::Trajectory t;
        t.append({100, -3, -4});
        t.append({105, 0, 0});
        CHECK(approx(t.length(), 5.0, 1e-12));
        CHECK(approx(t.durationMs(), 5.0, 1e-12));
        CHECK(approx(t.peakSpeed(), 1.0, 1e-12));
        CHECK(approx(t.peakSpeedTimeMs(), 2.5, 1e-12));
        CHECK(t.straightness() == 1.0);
    }

    // A path that only waits has no speed, and a step with no time is ignored for speed.
    {
        mq::Trajectory t;
        t.append({0, 1, 1});
        t.append({8, 1, 1});
        t.append({8, 9, 9}); // broken input: no time passes
        CHECK(t.peakSpeed() == 0.0);
        CHECK(t.peakSpeedTimeMs() == 0.0);
    }
}

// The rules every planned move must keep. Returns a description of the first broken rule,
// or nullptr when all hold.
static const char* brokenRule(const mq::Trajectory& t, mq::Point from, mq::Point to, const mq::MotionConfig& cfg) {
    if (t.empty()) return "empty";
    if (t.front().point() != from || t.front().tMs != 0.0) return "first sample is not the start at t = 0";
    if (t.back().point() != to) return "last sample is not the target";
    const auto& s = t.samples();
    for (size_t i = 1; i < s.size(); ++i) {
        if (s[i].tMs <= s[i - 1].tMs) return "timestamps do not strictly increase";
        if (s[i].point() == s[i - 1].point() && i + 1 != s.size()) return "duplicate pixel before the end";
    }
    if (t.durationMs() < mq::kMinDurationMs) return "shorter than the minimum duration";
    const double dt = 1000.0 / cfg.rateHz;
    if (t.durationMs() > mq::kMaxDurationMs + dt) return "longer than the maximum duration";

    // No teleports: no step is much longer than the fastest the cursor should ever go.
    // Minimum-jerk peaks at 1.875 × the average speed of the move.
    const double peak = 1.875 * mq::distance(mq::toVec2(from), mq::toVec2(to)) / t.durationMs();
    for (size_t i = 1; i < s.size(); ++i) {
        const double step = mq::distance(mq::toVec2(s[i - 1].point()), mq::toVec2(s[i].point()));
        if (step > peak * (s[i].tMs - s[i - 1].tMs) * 1.5 + 2.0) return "teleport";
    }
    return nullptr;
}

static void testFitts() {
    mq::MotionConfig cfg;
    CHECK(approx(mq::fittsIdBits(500, 16), std::log2(500.0 / 16 + 1), 1e-12));
    CHECK(mq::fittsIdBits(0, 16) == 0.0);

    // 80 + 130 · log2(500/16 + 1) ≈ 731.4 ms, and --speed divides it.
    const double t500 = 80 + 130 * std::log2(500.0 / 16 + 1);
    CHECK(approx(mq::fittsDurationMs(500, cfg), t500, 1e-9));
    cfg.speed = 2;
    CHECK(approx(mq::fittsDurationMs(500, cfg), t500 / 2, 1e-9));

    // Clamped at both ends.
    cfg.speed = 10;
    CHECK(mq::fittsDurationMs(10, cfg) == mq::kMinDurationMs);
    cfg.speed = 0.01;
    CHECK(mq::fittsDurationMs(3000, cfg) == mq::kMaxDurationMs);

    // A bigger target is quicker to hit.
    mq::MotionConfig big;
    big.targetSizePx = 64;
    CHECK(mq::fittsDurationMs(500, big) < mq::fittsDurationMs(500, mq::MotionConfig{}));

    // The minimum-jerk profile goes 0 → 1, is symmetric and never goes backwards.
    CHECK(mq::minimumJerk(0.0) == 0.0);
    CHECK(approx(mq::minimumJerk(1.0), 1.0, 1e-12));
    CHECK(approx(mq::minimumJerk(0.5), 0.5, 1e-12));
    CHECK(approx(mq::minimumJerk(0.2) + mq::minimumJerk(0.8), 1.0, 1e-12));
    bool monotonic = true;
    for (int i = 1; i <= 1000; ++i)
        if (mq::minimumJerk(i / 1000.0) < mq::minimumJerk((i - 1) / 1000.0)) monotonic = false;
    CHECK(monotonic);
}

static void testPlanMove() {
    const mq::MotionConfig cfg;

    // The rules hold for many seeds, distances and directions, also at negative coordinates.
    {
        const mq::Point starts[] = {{0, 0}, {-1500, 300}, {640, -200}};
        const double dirs[][2] = {{1, 0}, {0, 1}, {-1, 0}, {0, -1}, {0.6, 0.8}, {-0.8, 0.6}, {-0.6, -0.8}, {0.28, -0.96}};
        const int distances[] = {1, 5, 50, 500, 3000};
        const char* firstBroken = nullptr;
        int runs = 0;
        for (const mq::Point from : starts)
            for (const auto& d : dirs)
                for (const int dist : distances)
                    for (uint64_t seed = 0; seed < 40; ++seed) {
                        const mq::Point to{from.x + static_cast<int>(std::lround(d[0] * dist)),
                                           from.y + static_cast<int>(std::lround(d[1] * dist))};
                        mq::Rng rng(seed);
                        const char* broken = brokenRule(mq::planMove(from, to, cfg, rng), from, to, cfg);
                        if (broken && !firstBroken) {
                            firstBroken = broken;
                            std::fprintf(stderr, "  move (%d,%d) -> (%d,%d) seed %llu: %s\n", from.x, from.y,
                                         to.x, to.y, static_cast<unsigned long long>(seed), broken);
                        }
                        ++runs;
                    }
        CHECK(firstBroken == nullptr);
        CHECK(runs == 3 * 8 * 5 * 40);
    }

    // No distance, no move. One pixel still lands.
    {
        mq::Rng rng(1);
        CHECK(mq::planMove({7, 7}, {7, 7}, cfg, rng).empty());
        const mq::Trajectory one = mq::planMove({7, 7}, {8, 7}, cfg, rng);
        CHECK(brokenRule(one, {7, 7}, {8, 7}, cfg) == nullptr);
    }

    // Same seed, same path. Another seed, another path.
    {
        mq::Rng a(9), b(9), c(10);
        const auto ta = mq::planMove({0, 0}, {800, 450}, cfg, a);
        const auto tb = mq::planMove({0, 0}, {800, 450}, cfg, b);
        const auto tc = mq::planMove({0, 0}, {800, 450}, cfg, c);
        bool same = ta.size() == tb.size();
        for (size_t i = 0; same && i < ta.size(); ++i)
            same = ta.samples()[i].tMs == tb.samples()[i].tMs && ta.samples()[i].point() == tb.samples()[i].point();
        CHECK(same);
        CHECK(ta.durationMs() != tc.durationMs());
    }

    // Without jitter, longer moves never take less time (Fitts).
    {
        mq::MotionConfig steady = cfg;
        steady.jitterPct = 0;
        double prev = 0;
        bool grows = true;
        for (int dist = 1; dist <= 3000; dist += 7) {
            mq::Rng rng(1);
            const double d = mq::planMove({0, 0}, {dist, 0}, steady, rng).durationMs();
            if (d < prev) grows = false;
            prev = d;
        }
        CHECK(grows);
        mq::Rng r1(1), r2(1);
        CHECK(mq::planMove({0, 0}, {50, 0}, steady, r1).durationMs() < mq::planMove({0, 0}, {1000, 0}, steady, r2).durationMs());
    }

    // With jitter the average duration stays at Fitts's value.
    {
        double sum = 0;
        const int n = 2000;
        for (uint64_t seed = 0; seed < n; ++seed) {
            mq::Rng rng(seed);
            sum += mq::planMove({0, 0}, {500, 0}, cfg, rng).durationMs();
        }
        const double expected = mq::fittsDurationMs(500, cfg);
        CHECK(sum / n > expected * 0.97 && sum / n < expected * 1.03 + 8); // + up to one tick
    }

    // --duration fixes the length of the move (rounded up to the next tick).
    {
        mq::MotionConfig fixed = cfg;
        fixed.durationMs = 700;
        mq::Rng rng(3);
        const double d = mq::planMove({0, 0}, {300, 200}, fixed, rng).durationMs();
        CHECK(d >= 700 && d < 708);
    }

    // Samples sit on the update-rate grid, and the speed peaks in the middle (no skew yet).
    {
        mq::MotionConfig steady = cfg;
        steady.jitterPct = 0;
        steady.rateHz = 250;
        mq::Rng rng(4);
        const auto t = mq::planMove({0, 0}, {1000, 0}, steady, rng);
        bool onGrid = true;
        for (const auto& s : t.samples())
            if (!approx(s.tMs / 4.0, std::round(s.tMs / 4.0), 1e-9)) onGrid = false;
        CHECK(onGrid);
        const double at = t.peakSpeedTimeMs() / t.durationMs();
        CHECK(at > 0.4 && at < 0.6);
    }
}

int main() {
    testVersion();
    testVectors();
    testRect();
    testBezier();
    testArcLength();
    testRandom();
    testTrajectory();
    testFitts();
    testPlanMove();
    if (failures == 0) std::puts("all tests passed");
    return failures == 0 ? 0 : 1;
}
