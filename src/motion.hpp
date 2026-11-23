#pragma once

#include "geometry.hpp"
#include "random.hpp"
#include "trajectory.hpp"

namespace mq {

// Everything that shapes one move. All values are final: profiles and command-line flags are
// merged into this before planning.
struct MotionConfig {
    double durationMs = 0.0;    // fixed duration of the main movement; 0 = automatic (Fitts)
    double speed = 1.0;         // multiplier on the automatic duration (2 = twice as fast)
    double jitterPct = 15.0;    // random variation of the automatic duration (log-normal)
    double targetSizePx = 16.0; // target width W for Fitts's law
    double rateHz = 125.0;      // cursor update rate

    double curvature = 0.30;      // how much paths bend, 0 = straight
    double sCurveChancePct = 5.0; // chance that the path bends one way, then the other

    // Fitts's law constants: T = a + b · log2(D / W + 1).
    double fittsAMs = 80.0;
    double fittsBMs = 130.0;
};

constexpr double kMinDurationMs = 60.0;
constexpr double kMaxDurationMs = 5000.0;

// Fitts's index of difficulty in bits: log2(D / W + 1).
double fittsIdBits(double distancePx, double targetSizePx);

// Automatic duration of a move of `distancePx` without jitter: Fitts's law divided by
// cfg.speed, clamped to [kMinDurationMs, kMaxDurationMs].
double fittsDurationMs(double distancePx, const MotionConfig& cfg);

// Minimum-jerk position profile: 0 at τ = 0, 1 at τ = 1, zero speed and acceleration at both
// ends, bell-shaped speed in between (peak 1.875 at τ = 0.5).
inline double minimumJerk(double tau) {
    const double t3 = tau * tau * tau;
    return t3 * (10.0 + tau * (-15.0 + 6.0 * tau));
}

// Distance used to scale the bend of a path. Equal to the distance up to 300 px, then it
// grows more slowly, so long moves are not absurdly curved.
double bendDistancePx(double distancePx);

// The shape of a move from `a` to `b`: a cubic Bézier whose control points lie at 30 % and
// 70 % of the way, pushed sideways by a random amount scaled by cfg.curvature. Both are
// pushed to the same side (one smooth arc), except with probability cfg.sCurveChancePct,
// when the second one goes to the other side (a gentle S).
CubicBezier planPath(Vec2 a, Vec2 b, const MotionConfig& cfg, Rng& rng);

// Plans the cursor path from `from` to `to`.
//
// The cursor follows planPath() with the minimum-jerk speed profile. Positions are taken by
// arc length, so the speed along a curve is right. The first sample is `from` at t = 0, the last one is exactly `to`. Samples are taken every
// 1000 / rateHz ms. A sample on the same pixel as the previous one is dropped (time still
// advances), except the final one, which always marks the end of the move. A move shorter
// than 1 px returns an empty trajectory.
Trajectory planMove(Point from, Point to, const MotionConfig& cfg, Rng& rng);

} // namespace mq
