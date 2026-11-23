#include "trajectory.hpp"

namespace mq {

double Trajectory::durationMs() const {
    return samples_.size() < 2 ? 0.0 : samples_.back().tMs - samples_.front().tMs;
}

double Trajectory::length() const {
    double total = 0.0;
    for (size_t i = 1; i < samples_.size(); ++i)
        total += distance(toVec2(samples_[i - 1].point()), toVec2(samples_[i].point()));
    return total;
}

double Trajectory::straightness() const {
    const double len = length();
    if (len <= 0.0) return 1.0;
    return distance(toVec2(samples_.front().point()), toVec2(samples_.back().point())) / len;
}

double Trajectory::stepSpeed(size_t i) const {
    const Sample& a = samples_[i];
    const Sample& b = samples_[i + 1];
    const double dt = b.tMs - a.tMs;
    if (dt <= 0.0) return 0.0;
    return distance(toVec2(a.point()), toVec2(b.point())) / dt;
}

size_t Trajectory::fastestStep() const {
    size_t best = samples_.size();
    double bestSpeed = 0.0;
    for (size_t i = 0; i + 1 < samples_.size(); ++i) {
        const double v = stepSpeed(i);
        if (v > bestSpeed) {
            bestSpeed = v;
            best = i;
        }
    }
    return best;
}

double Trajectory::peakSpeed() const {
    const size_t i = fastestStep();
    return i < samples_.size() ? stepSpeed(i) : 0.0;
}

double Trajectory::peakSpeedTimeMs() const {
    const size_t i = fastestStep();
    if (i >= samples_.size()) return 0.0;
    return (samples_[i].tMs + samples_[i + 1].tMs) / 2.0 - samples_.front().tMs;
}

} // namespace mq
