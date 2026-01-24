# CLAUDE.md

Guidance for Claude Code working with this repository.

## Project Overview

**Conformality** is a C++20 conformal mapping toolset implementing the Fornberg MC method for multiply-connected domains. Uses pure composition with a single orchestrator pattern.

## Prerequisites

**Required (macOS):**
```bash
# MacPorts
sudo port install fftw-3 eigen3 glfw

# Or Homebrew
brew install fftw eigen glfw
```

**Verify installation:**
```bash
ls /opt/local/lib/libfftw3.* || ls /usr/local/lib/libfftw3.*   # FFTW
ls /opt/local/include/eigen3 || ls /usr/local/include/eigen3   # Eigen
```

## Quick Commands

Use skills for common tasks:
- `/build` - Build the project (add "clean" or "reconfigure" as needed)
- `/test` - Run tests (optionally with gtest filter)
- `/gui` - Launch the GUI application
- `/project-status` - Check GitHub issues and milestones

### Manual Build
```bash
cmake -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build -j$(sysctl -n hw.ncpu)
ctest --test-dir build         # Run tests
./build/conformality_gui       # Run GUI
```

## Directory Structure

```
src/core/      - ConformalMap, Types, StatusManager
src/methods/   - TheodorsenMethod, FornbergMC, PMatrixBuilder
src/domains/   - Domain, Boundary, FornbergCanonicalDomain
src/numerics/  - FFTWWrapper, Grid, RootFinder, CGSolver, Polyval
src/gui/       - ImGui-based GUI components
tests/         - Google Test suite
external/      - Third-party (imgui, implot, spdlog)
design/        - Design docs and MATLAB reference
```

## Key Entry Points

| To understand... | Start here |
|------------------|------------|
| Overall architecture | `src/core/ConformalMap.h` - single orchestrator |
| Current algorithm work | `src/methods/FornbergMC.h` + `PMatrixBuilder.h` |
| Algorithm configuration | `src/methods/FornbergMCConfiguration.h` |
| Domain representation | `src/domains/Domain.h` → `MultiplyConnectedDomain` |
| Numerical core | `src/numerics/FFTWWrapper.h`, `CGSolver.h` |
| GUI application | `src/gui/Application.cpp` → `MainWindow.cpp` |

## C++ Style

**C++20** - Guidelines, not rules. Break them if it improves readability.

| Element | Convention | Example |
|---------|------------|---------|
| Variables | `snake_case` | `point_count` |
| Classes | `PascalCase` | `ConformalMap` |
| Functions | `camelCase` | `computeCoefficients` |
| Private members | `m_` prefix | `m_value` |
| Pointers | `p_` prefix | `mp_data` (member ptr) |

**Formatting:** 120 char lines, braces on own lines, curly bracket init preferred

## Key Patterns

- **Error handling:** `std::invalid_argument` (validation), `std::runtime_error` (runtime)
- **Ownership:** `shared_ptr` for domains/methods, `unique_ptr` for internals
- **Libraries:** Eigen (linear algebra), FFTW via `FFTWWrapper`, spdlog (logging)
- **Test builds:** `TESTING` preprocessor define enabled for test executable

## Reference Documentation

Detailed information in `.claude/references/conformality/`:

| Topic | File |
|-------|------|
| Fornberg algorithm | `fornberg-method.md` |
| Theodorsen method (on hold) | `theodorsen-method.md` |
| FFTW & numerics | `fftw-numerics.md` |
| Testing guidelines | `testing.md` |
| Architecture | `architecture.md` |
| Project phases | `project-phases.md` |

## Gotchas

- **FFTW wisdom:** `fftw.wisdom` file generated in working directory at runtime
- **ImGui config:** `imgui.ini` file generated in working directory at runtime
- **macOS Metal:** GUI uses Metal backend; won't build on Linux/Windows without porting
- **Disabled tests:** Some tests prefixed `DISABLED_` until implementation complete
- **MATLAB reference:** `design/fornberg/fornmc/` is a symlink to reference implementation
- **Design directory:** Has its own `.git` for development (will be removed before publish)

## Troubleshooting

| Problem | Solution |
|---------|----------|
| `Unable to find fftw3.h` | Install FFTW: `port install fftw-3` or `brew install fftw` |
| `Unable to find Eigen3` | Install Eigen: `port install eigen3` or `brew install eigen` |
| `GLFW not found` | Install GLFW: `port install glfw` or `brew install glfw` |
| CMake can't find libs | Check paths in CMakeLists.txt match your install (`/opt/local` vs `/usr/local`) |
| Tests crash on FFT | Delete `fftw.wisdom` and re-run (wisdom file may be stale) |

## Notes

- Strip trailing whitespace from all files
- License: AGPL v3
