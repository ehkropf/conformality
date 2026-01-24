# CLAUDE.md

Guidance for Claude Code working with this repository.

# Claude Project Guide

**Conformality** is a C++20 conformal mapping toolset implementing the Fornberg MC method for multiply-connected domains. Uses pure composition with a single orchestrator pattern.

**Current Phase:** Phase 1 - Theodorsen method implementation with GUI
**Project Overview:** `design/high_level_overview.md`
**Phase 1 Documentation:** `design/phase1/` directory
**Map Design Discussion:** `design/phase1/map_hierarchy_design_discussion.md` (led to major refactor)

### Current Implementation Status
- ✅ Core domain and mapping architecture completed
- ✅ Theodorsen method for simply connected domains
- ✅ ImGui-based GUI with Metal rendering (macOS)
- ✅ Comprehensive test suite with Google Test
- ✅ FFTW integration for Laurent series computation
- 🔄 Real-time parameter adjustment and visualization
- 📋 Planned: Numerical tolerance management refactor

## Project Architecture

**Conformality** is a C++20 conformal mapping toolset using pure composition with a single orchestrator pattern.

### Core Components

#### 1. Conformal Mapping System
- **`ConformalMap`** - Single orchestrator class coordinating domain and method
- **`ConformalMapMethod`** - Abstract base for mapping algorithms
- **`TheodorsenMethod`** - Implements Theodorsen's method for simply connected domains

#### 2. Domain System
- **`Domain`** - Abstract base class for geometric domains (elliptical, circular, starlike)
- **`Boundary`** - Manages domain boundaries with parameterized components
- **`BoundaryComponent`** - Individual boundary curve segments with sampling

#### 3. Numerical Infrastructure
- **`FFTWWrapper`** - FFTW library integration for Laurent series coefficients
- **`Grid`** - Grid generation and conformal grid visualization
- **`RootFinder`** - Numerical root finding with configurable tolerances
- **`Types`** - Core mathematical types and complex number utilities

#### 4. GUI System (ImGui + Metal)
- **`Application`** - Main application controller and event loop
- **`MainWindow`** - Primary GUI window with dual-domain visualization
- **`VisualizationPanel`** - Real-time conformal mapping display
- **`MetalRenderer`** - Metal-based rendering backend for macOS
- **`GuiController`** - GUI state management and parameter controls

### Directory Structure
- `src/` - Core implementation (domains, methods, numerical tools)
- `src/gui/` - ImGui-based GUI components
- `tests/` - Google Test suite with comprehensive coverage
- `external/` - Third-party libraries (imgui, implot)
- `design/` - Phase documentation and architectural decisions

## Dependencies & Libraries

- **Third-party included libraries**
  - imgui
  - implot
  - spdlog
- **Libraries that must be pre-installed**
  - FFTW3 (Fast Fourier Transform library, needs compiled for target system)
  - GLFW (windowing library)
  - Metal/MetalKit (macOS graphics - automatically available on macOS)
  - Cocoa/IOKit/CoreVideo/QuartzCore (macOS frameworks)

## Common patterns & Utilities

- **Error handling**
  - It's good to throw exceptions on errors (except when it isn't)
  - Errors should have good messages
- **Logging**
  - StatusManager will handle collecting messages and logging
  - The `spdlog` package will be used for logging

## FFTW & Numerical Methods

### FFTW Coefficient Ordering
- **Forward FFT output order:** `[0, 1, 2, ..., N/2-1, -N/2, -N/2+1, ..., -1]`
- **Normalization:** FFTW does not normalize by default - divide by N for proper DFT coefficients
- **Laurent Series Mapping:** For conformal mapping applications:
  - Positive frequencies (indices 1 to N/2-1) → z^k terms
  - Negative frequencies (indices N/2 to N-1) → z^{-j} terms where j = N-k
  - Index 0 → constant term
  - Index N-1 → z^{-1} coefficient (dominant for conformal maps)

### Theodorsen Method Specifics
- **Series form:** f(z) = a₀ + a₋₁z + a₁/z + a₂/z² + ...
- **Numerical stability:** High-order positive powers z^k can cause explosion for |z| > 1
- **Practical evaluation:** Use only constant, linear (z), and first few inverse power (1/z^k) terms
- **Convergence issues:** Extreme ellipse aspect ratios (>2:1) may not converge within default iteration limits

## Build & Development Commands

### Building the Project
```bash
# MacPorts
sudo port install fftw-3 eigen3 glfw

# Or Homebrew
brew install fftw eigen glfw
```

**Verify installation:**
```bash
# From build directory
./test_exec
```

## Directory Structure

### Key Build Targets
- Main executable: `conformality` (GUI application)
- Test executable: `test_exec` (Google Test suite)

## Documentation & Diagrams
- Automatic class diagrams are generated with `clang-uml`
- Configuration file is `.clang-uml` in the project root directory
- Auto generated diagrams are available in docs/diagrams. Claude may generate any diagrams in the config file not already in this directory.

## Project Management
- GitHub integration available via `gh` tool
- Can use milestones and issues for planning discussions

## C++ Style Guide

This project uses c++20.

*Note: These are guidelines and can be broken if it improves readability.*

### General Formatting
- **Line length:** Keep lines under 120 characters
- **Braces:** Always on their own lines (exception: empty functions can use `{}` on one line)
- **Spacing:** Use spaces around operators (`+`, `-`) and after commas

### Naming Conventions
- **Variables:** `snake_case`
- **Classes:** `PascalCase`
- **Functions:** `camelCase`
- **Private members:** Prefix with `m_`
- **Pointers:** Prefix with `p_` (member pointers: `mp_`)

### Function Formatting
- **Arguments:** All on same line unless exceeding line limit
- **Continuation:** Align under first character of previous argument, or indent 8 spaces if no previous argument
- **Closing parenthesis:** Same line as last argument

### Constructor Initialization
- **Colon:** On separate line to start the list
- **Members:** Each on new line starting with `,` aligned under the colon

### Object Initialization
- **Preferred:** Curly brackets `{}`
- **Alternative:** Parentheses or assignment when it improves readability

### Example
```cpp
class MyClass
{
private:
    int m_value;
    std::unique_ptr<int> mp_data;

public:
    MyClass(int value, std::unique_ptr<int> data)
        : m_value{value}
        , mp_data{std::move(data)}
    {
    }

    void myFunction(int param1, int param2)
    {
        const int result = param1 + param2;
        // Function body
    }
};
```

## Key Entry Points

### Memory & Resources

- Use smart pointers (`std::unique_ptr`, `std::shared_ptr`) for dynamic memory
- Prefer `std::vector` over raw arrays (unless compile-time size is known)
- Avoid global variables

### Code Quality

- Always include necessary headers
- Use `const` and `constexpr` extensively (declare const-ness unless value will change)
- Prefer RAII for resource management
- When possible, always go for the best big-O number for performance

| Element | Convention | Example |
|---------|------------|---------|
| Variables | `snake_case` | `point_count` |
| Classes | `PascalCase` | `ConformalMap` |
| Functions | `camelCase` | `computeCoefficients` |
| Private members | `m_` prefix | `m_value` |
| Pointers | `p_` prefix | `mp_data` (member ptr) |

**Formatting:** 120 char lines, braces on own lines, curly bracket init preferred

### Test Coverage Areas
- **Boundary System:** Component creation, parameterization, sampling
- **Domain Implementations:** Elliptical, circular, and starlike domains
- **Conformal Mapping:** Theodorsen method with internal/external mappings
- **Numerical Methods:** FFT wrapper, root finding, coefficient computation
- **Status Management:** Error handling and convergence monitoring

### Test Organization
- Tests are comprehensive and cover all major system components
- Each test file focuses on a specific module or class
- Integration tests verify component interactions
- Numerical tests include tolerance and convergence validation
