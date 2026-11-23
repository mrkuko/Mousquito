# Mousquito

[![CI](https://github.com/mrkuko/Mousquito/actions/workflows/ci.yml/badge.svg)](https://github.com/mrkuko/Mousquito/actions/workflows/ci.yml)
[![License: GPL v2](https://img.shields.io/badge/License-GPL_v2-blue.svg)](LICENSE)
![Platforms](https://img.shields.io/badge/platform-Windows%20%7C%20Linux%20%7C%20macOS-lightgrey)
![C++17](https://img.shields.io/badge/C%2B%2B-17-00599C)

Mousquito moves the mouse cursor the way a person would: curved paths, a natural speed curve, a slight overshoot and a small correction. It is the sister tool of [Pasteguise](https://github.com/mrkuko/Pasteguis): Pasteguise does human-like typing, Mousquito does human-like pointing.

Useful for screencasts and product demos, UI test automation, accessibility macros, presentations, and keeping a machine awake.

> **Status: early development.** The project builds and its tests run, but it does not move the mouse yet. The options below are the planned interface.

## Planned features

- Cross-platform: Windows, Linux (X11) and macOS
- Curved paths with a bell-shaped speed curve, not straight robotic lines
- Duration from Fitts's law: longer and smaller targets take longer
- Overshoot, undershoot and small corrective moves, always landing exactly on the target
- Subtle correlated hand tremor and occasional mid-move slow-downs ("hiccups")
- Profiles: `precise`, `casual`, `hurried`, `slow`, `tired`
- Reproducible runs with `--seed`
- `--dry-run` statistics and `--export` of the planned path to CSV or SVG
- `--buzz`: hover restlessly around a spot, like a mosquito
- Hold **Esc** to swat it. Move the mouse yourself and it quits at once.

## Build

Requires CMake 3.16+ and a C++17 compiler.

```sh
cmake -S . -B build
cmake --build build --config Release
ctest --test-dir build -C Release
```

| Platform | Extra requirements |
|----------|--------------------|
| Windows  | Visual Studio 2019+ or MinGW |
| Linux    | `libx11-dev libxtst-dev` (once the X11 backend lands) |
| macOS    | Xcode command line tools, and **Accessibility** permission for your terminal |

## Planned usage

```sh
mousquito 800,450                         # glide to 800,450 with the default "casual" profile
mousquito 100,100 900,300 400,700         # visit several points in order
mousquito --by -200,100                   # relative move
mousquito 1200,700 --profile hurried      # faster, sloppier, more overshoot
mousquito 800,450 --dry-run --seed 7      # print the plan, move nothing
mousquito 800,450 --seed 7 --export path.svg
mousquito --buzz 60                       # hover around the cursor for 60 s
```

Coordinates are global desktop pixels (`x,y`). Negative values are allowed on multi-monitor setups.

## What it is not

Mousquito is not a tool for bypassing CAPTCHAs, bot detection or anti-cheat systems, and it makes no "undetectable" claims. The profiles are parameter presets, not imitations of any real person's movement.

## Platform notes

- **Linux:** X11/XWayland only (XTest). Pure Wayland sessions are not supported.
- **Windows:** input cannot be sent to elevated (administrator) windows from a normal process.
- **macOS:** needs Accessibility permission.

## License

Copyright (C) 2026 Emanuel Cisar

Mousquito is free software, licensed under the [GNU General Public License v2.0](LICENSE). You may use, modify and redistribute it, but you must keep this copyright notice and license, and distribute derivative works under the same license with source available.
