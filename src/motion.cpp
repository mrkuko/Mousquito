#include "motion.hpp"

#include <algorithm>
#include <cmath>

namespace mq {

double fittsIdBits(double distancePx, double targetSizePx) {
    const double w = std::max(targetSizePx, 1.0);
    return std::log2(std::max(distancePx, 0.0) / w + 1.0);
}

double fittsDurationMs(double distancePx, const MotionConfig& cfg) {
    const double speed = cfg.speed > 0.0 ? cfg.speed : 1.0;
    const double t = (cfg.fittsAMs + cfg.fittsBMs * fittsIdBits(distancePx, cfg.targetSizePx)) / speed;
    return std::clamp(t, kMinDurationMs, kMaxDurationMs);
}

namespace {

// Duration of the main movement: --duration if given, otherwise Fitts × jitter, clamped.
// Always draws the jitter, so a fixed duration does not shift later random decisions.
double mainDurationMs(double distancePx, const MotionConfig& cfg, Rng& rng) {
    const double jitter = rng.logNormalMean1(cfg.jitterPct / 100.0);
    if (cfg.durationMs > 0.0) return cfg.durationMs;
    return std::clamp(fittsDurationMs(distancePx, cfg) * jitter, kMinDurationMs, kMaxDurationMs);
}

} // namespace

Trajectory planMove(Point from, Point to, const MotionConfig& cfg, Rng& rng) {
    Trajectory traj;
    const Vec2 a = toVec2(from);
    const Vec2 b = toVec2(to);
    const double dist = distance(a, b);
    if (dist < 1.0) return traj;

    const double durationMs = mainDurationMs(dist, cfg, rng);
    const double dtMs = 1000.0 / (cfg.rateHz > 0.0 ? cfg.rateHz : 125.0);

    // Samples on the fixed grid t_k = k · dt. The last one is the first grid tick at or after
    // the end of the move, and it is exactly the target.
    const long lastTick = std::max(1L, static_cast<long>(std::ceil(durationMs / dtMs - 1e-9)));
    traj.append({0.0, from.x, from.y});
    Point prev = from;
    for (long k = 1; k < lastTick; ++k) {
        const double tMs = static_cast<double>(k) * dtMs;
        const double s = minimumJerk(tMs / durationMs);
        const Point p = roundToPoint(a + (b - a) * s);
        if (p == prev) continue;
        traj.append({tMs, p.x, p.y});
        prev = p;
    }
    traj.append({static_cast<double>(lastTick) * dtMs, to.x, to.y});
    return traj;
}

} // namespace mq
