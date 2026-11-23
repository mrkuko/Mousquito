#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

namespace mq {

// A point or direction in continuous desktop space (pixels, may be fractional).
struct Vec2 {
    double x = 0.0;
    double y = 0.0;

    Vec2() = default;
    Vec2(double x_, double y_) : x(x_), y(y_) {}

    Vec2 operator+(Vec2 o) const { return {x + o.x, y + o.y}; }
    Vec2 operator-(Vec2 o) const { return {x - o.x, y - o.y}; }
    Vec2 operator*(double k) const { return {x * k, y * k}; }
    Vec2 operator/(double k) const { return {x / k, y / k}; }
    Vec2& operator+=(Vec2 o) { x += o.x; y += o.y; return *this; }
    Vec2& operator-=(Vec2 o) { x -= o.x; y -= o.y; return *this; }

    double length() const { return std::hypot(x, y); }
};

inline Vec2 operator*(double k, Vec2 v) { return v * k; }
inline double dot(Vec2 a, Vec2 b) { return a.x * b.x + a.y * b.y; }
inline double distance(Vec2 a, Vec2 b) { return (b - a).length(); }

// Unit vector in the direction of `v`, or (0,0) when `v` has no length.
inline Vec2 normalized(Vec2 v) {
    const double len = v.length();
    return len > 0.0 ? v / len : Vec2{};
}

// Unit vector perpendicular to `dir` (rotated +90°), or (0,0) when `dir` has no length.
inline Vec2 unitNormal(Vec2 dir) {
    const Vec2 n = normalized(dir);
    return {-n.y, n.x};
}

// A whole desktop pixel.
struct Point {
    int x = 0;
    int y = 0;

    bool operator==(Point o) const { return x == o.x && y == o.y; }
    bool operator!=(Point o) const { return !(*this == o); }
};

inline Vec2 toVec2(Point p) { return {static_cast<double>(p.x), static_cast<double>(p.y)}; }

// Nearest whole pixel.
inline Point roundToPoint(Vec2 v) {
    return {static_cast<int>(std::lround(v.x)), static_cast<int>(std::lround(v.y))};
}

// A pixel rectangle. May start at negative coordinates (monitors left of or above the primary).
struct Rect {
    int x = 0;
    int y = 0;
    int w = 0;
    int h = 0;

    bool empty() const { return w <= 0 || h <= 0; }

    bool contains(Point p) const {
        return !empty() && p.x >= x && p.x < x + w && p.y >= y && p.y < y + h;
    }

    // Nearest pixel inside the rectangle. An empty rectangle returns `p` unchanged.
    Point clamp(Point p) const {
        if (empty()) return p;
        return {std::clamp(p.x, x, x + w - 1), std::clamp(p.y, y, y + h - 1)};
    }
};

// Cubic Bézier curve P0 → P3 with control points P1 and P2.
struct CubicBezier {
    Vec2 p0, p1, p2, p3;

    Vec2 at(double t) const {
        const double u = 1.0 - t;
        return p0 * (u * u * u) + p1 * (3.0 * u * u * t) + p2 * (3.0 * u * t * t) + p3 * (t * t * t);
    }

    // First derivative dB/dt.
    Vec2 derivative(double t) const {
        const double u = 1.0 - t;
        return (p1 - p0) * (3.0 * u * u) + (p2 - p1) * (6.0 * u * t) + (p3 - p2) * (3.0 * t * t);
    }
};

// Maps arc length s (pixels travelled along the curve) to the Bézier parameter t and back.
//
// Bézier t is not uniform in space: equal steps in t give unequal steps along the path.
// Moving at a given speed therefore needs positions by arc length. The table samples the
// curve at `samples` evenly spaced t values, sums the chord lengths, and inverts by
// binary search with linear interpolation between neighbouring samples.
class ArcLengthTable {
public:
    static constexpr int kDefaultSamples = 256;
    static constexpr double kZeroLength = 1e-9; // pixels

    explicit ArcLengthTable(const CubicBezier& curve, int samples = kDefaultSamples) {
        const int n = std::max(samples, 2);
        cumulative_.resize(static_cast<size_t>(n));
        cumulative_[0] = 0.0;
        Vec2 prev = curve.at(0.0);
        for (int i = 1; i < n; ++i) {
            const Vec2 cur = curve.at(static_cast<double>(i) / (n - 1));
            cumulative_[static_cast<size_t>(i)] = cumulative_[static_cast<size_t>(i - 1)] + distance(prev, cur);
            prev = cur;
        }
        // A curve that does not move still sums to ~1e-15 from rounding in at(). Call it 0.
        if (cumulative_.back() < kZeroLength) std::fill(cumulative_.begin(), cumulative_.end(), 0.0);
    }

    // Total length of the curve in pixels.
    double length() const { return cumulative_.back(); }

    // Arc length from the start of the curve to parameter t (clamped to [0, 1]).
    double lengthAtT(double t) const {
        const double pos = std::clamp(t, 0.0, 1.0) * lastIndex();
        const size_t i = std::min(static_cast<size_t>(pos), cumulative_.size() - 2);
        const double frac = pos - static_cast<double>(i);
        return cumulative_[i] + (cumulative_[i + 1] - cumulative_[i]) * frac;
    }

    // Parameter t at arc length s (clamped to [0, length()]). A zero-length curve returns 0.
    double tAtLength(double s) const {
        const double total = length();
        if (total <= 0.0 || s <= 0.0) return 0.0;
        if (s >= total) return 1.0;

        // First sample with cumulative length >= s. It is at index >= 1 because s > 0.
        const auto it = std::lower_bound(cumulative_.begin(), cumulative_.end(), s);
        const size_t hi = static_cast<size_t>(it - cumulative_.begin());
        const size_t lo = hi - 1;
        const double segment = cumulative_[hi] - cumulative_[lo];
        const double frac = segment > 0.0 ? (s - cumulative_[lo]) / segment : 0.0;
        return (static_cast<double>(lo) + frac) / lastIndex();
    }

private:
    double lastIndex() const { return static_cast<double>(cumulative_.size() - 1); }

    std::vector<double> cumulative_;
};

} // namespace mq
