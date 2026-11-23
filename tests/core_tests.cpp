#include <cmath>
#include <cstdio>
#include <cstring>

#include "geometry.hpp"
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

int main() {
    testVersion();
    testVectors();
    testRect();
    testBezier();
    testArcLength();
    testRandom();
    testTrajectory();
    if (failures == 0) std::puts("all tests passed");
    return failures == 0 ? 0 : 1;
}
