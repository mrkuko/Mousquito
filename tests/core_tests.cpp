#include <cmath>
#include <cstdio>
#include <cstring>

#include "geometry.hpp"
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

int main() {
    testVersion();
    testVectors();
    testRect();
    testBezier();
    testArcLength();
    if (failures == 0) std::puts("all tests passed");
    return failures == 0 ? 0 : 1;
}
