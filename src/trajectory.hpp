#pragma once

#include <cstddef>
#include <vector>

#include "geometry.hpp"

namespace mq {

// One cursor position at a moment of the move. Time is relative to the start of the move.
struct Sample {
    double tMs = 0.0;
    int x = 0;
    int y = 0;

    Point point() const { return {x, y}; }
};

// The planned path of the cursor: samples in time order, ready to be played or inspected.
class Trajectory {
public:
    // Adds a sample at the end. The caller keeps timestamps strictly increasing.
    void append(Sample s) { samples_.push_back(s); }

    const std::vector<Sample>& samples() const { return samples_; }
    bool empty() const { return samples_.empty(); }
    size_t size() const { return samples_.size(); }
    const Sample& front() const { return samples_.front(); }
    const Sample& back() const { return samples_.back(); }

    // Time from the first to the last sample. 0 with fewer than two samples.
    double durationMs() const;

    // Distance travelled along the samples, in pixels.
    double length() const;

    // Straight-line distance from start to end divided by length(): 1 for a straight path,
    // lower for curved ones. 1 for a path that does not move.
    double straightness() const;

    // Highest speed between two consecutive samples, in pixels per millisecond.
    double peakSpeed() const;

    // When that speed happens: the midpoint in time of the fastest step. 0 if there is none.
    double peakSpeedTimeMs() const;

private:
    // Index i of the fastest step samples_[i] -> samples_[i + 1], or size() if there is none.
    size_t fastestStep() const;
    double stepSpeed(size_t i) const;

    std::vector<Sample> samples_;
};

} // namespace mq
