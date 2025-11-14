# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

# Conformality Project Guide

## Project Status & Documentation

**Current Branch:** `fornmc`
**Current Phase:** Phase 1 - Fornberg MC Implementation
**Project Overview:** `design/overview.md`
**Roadmap:** `design/roadmap.md`
**Fornberg Documentation:** `design/fornberg/` directory

### Current Implementation Status

**Active Development - Fornberg MC Method:**
- ✅ Core infrastructure complete (FornbergMC, PMatrixBuilder, CGSolver, FornbergCanonicalDomain)
- ✅ Configuration system with ~30 tunable parameters
- ✅ Test framework in place
- 🔄 Core algorithm methods ~60% complete (~36 TODO/FIXME markers)
- 🎯 **Current Goal:** Complete core implementation and validate using constructive tests

**On Hold - Theodorsen Method:**
- ⚠️ Formulation error requiring additional research
- ✅ GUI and infrastructure exist
- 📋 Status: On hold until Fornberg work complete

**Completed - Foundation:**
- ✅ Core domain and mapping architecture
- ✅ ImGui-based GUI with Metal rendering (macOS)
- ✅ Comprehensive test suite with Google Test
- ✅ FFTW integration for FFT operations
- ✅ Numerical infrastructure (Grid, RootFinder, FFTWWrapper)

## Project Architecture

**Conformality** is a C++20 conformal mapping toolset using pure composition with a single orchestrator pattern.

### Design Philosophy

**Pure Composition Over Inheritance**
- Single orchestrator pattern (`ConformalMap`)
- Clear ownership and lifetime management
- Easy to reason about data flow
- Flexible for adding new methods

### Core Components

#### 1. Conformal Mapping System
- **`ConformalMap`** - Single orchestrator coordinating domain and method
- **`ConformalMapMethod`** - Abstract base for mapping algorithms

**Implemented Methods:**
- **`TheodorsenMethod`** - Simply-connected domains (ON HOLD - needs research)
- **`FornbergMC`** - Multiply-connected domains (IN PROGRESS)

#### 2. Domain System
- **`Domain`** - Abstract base class for geometric domains
- **`MultiplyConnectedDomain`** - Extends Domain for m-connected regions
- **`FornbergCanonicalDomain`** - Canonical domain (unit disk with circular holes)
- **`Boundary`** - Manages domain boundaries with parameterized components
- **`BoundaryComponent`** - Individual boundary curve segments with sampling

#### 3. Numerical Infrastructure
- **`FFTWWrapper`** - FFTW library integration for FFT operations
- **`CGSolver`** - Custom conjugate gradient solver (ported from MATLAB)
- **`PMatrixBuilder`** - P matrix construction for analyticity conditions
- **`Polyval`** - Polynomial evaluation utilities
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
- `src/core/` - ConformalMap, Types, StatusManager
- `src/methods/` - TheodorsenMethod, FornbergMC, PMatrixBuilder
- `src/domains/` - Domain, Boundary, FornbergCanonicalDomain
- `src/numerics/` - FFTWWrapper, Grid, RootFinder, CGSolver, Polyval
- `src/gui/` - ImGui-based GUI components
- `tests/` - Google Test suite with comprehensive coverage
- `external/` - Third-party libraries (imgui, implot, spdlog)
- `design/` - Design documentation and experiments
  - `design/fornberg/` - Fornberg method documentation
  - `design/lacompare/` - Linear algebra benchmarking (Eigen vs Armadillo)
  - `design/cpp_flag_examples/` - Expression template experiments
  - `design/plotting/` - ImGui plotting validation

## Method Overview

### Fornberg MC Method (Active Development)

**Purpose:** Conformal mapping for multiply-connected domains (m ≥ 2)

**Algorithm:** Newton iteration with analyticity conditions
- Solves for boundary correspondences S_ν(θ) and conformal moduli (c_ν, ρ_ν)
- Uses Fourier series representation with FFT
- Custom CG solver for linear systems
- Special optimization for annulus case (m=2)

**Key Files:**
- `src/methods/FornbergMC.{h,cpp}` - Main algorithm (~22 TODO markers)
- `src/methods/FornbergMCConfiguration.h` - Comprehensive configuration
- `src/methods/PMatrixBuilder.{h,cpp}` - Matrix construction (~3 TODO markers)
- `src/numerics/CGSolver.{h,cpp}` - Conjugate gradient solver (~4 FIXME markers)
- `src/domains/FornbergCanonicalDomain.{h,cpp}` - Canonical domain (~2 FIXME markers)
- `tests/fornberg_mc.cpp` - Test suite (~2 TODO markers)
- `tests/cgsolver.cpp` - CG solver tests
- `tests/pmatrixbuilder.cpp` - Matrix builder tests

**MATLAB Reference:** `design/fornberg/fornmc/` symlink (points to ~/Source/gdrive/fornmc/matlab)

**Status:** See `design/fornberg/current_status.md` for detailed implementation status

**Critical Missing Implementations:**
1. `formSystem()` - System matrix and RHS construction
2. `solveSystem()` - CG solver integration
3. `newtonUpdate()` - Parameter update logic
4. `sampleBoundaries()` - Boundary sampling
5. `computeFourierCoefficients()` - FFT-based coefficient computation
6. `redistributeBoundaryParameters()` - Adaptive point redistribution
7. `inverseMap()` - Inverse mapping via Newton's method
8. Logging integration throughout (~16 TODO markers)

### Theodorsen Method (On Hold)

**Purpose:** Fast conformal mapping for simply-connected domains

**Status:** ⚠️ Has formulation error requiring additional research

**Algorithm:** FFT-based Laurent series computation

**Key Challenges:**
- Numerical stability for |z| > 1
- Extreme aspect ratios (>2:1) convergence issues
- Current implementation has unresolved formulation bug

**Files:**
- `src/methods/TheodorsenMethod.{h,cpp}`
- `tests/theodorsen_method_tests.cpp`

**Note:** Will return to this after Fornberg implementation complete

## Dependencies & Libraries

### Third-party Included Libraries
- **imgui** - Immediate mode GUI
- **implot** - Real-time plotting
- **spdlog** - Fast logging library

### Pre-installed System Dependencies
- **FFTW3** - Fast Fourier Transform (must be compiled for target system)
- **Eigen** - Linear algebra (header-only)
- **GLFW** - Windowing library
- **Google Test** - Testing framework

### Platform Dependencies (macOS)
- **Metal/MetalKit** - GPU rendering (automatically available)
- **Cocoa/IOKit/CoreVideo/QuartzCore** - System frameworks

## Common Patterns & Utilities

### Error Handling
- Throw exceptions for errors with descriptive messages
- Use `std::invalid_argument` for validation failures
- Use `std::runtime_error` for runtime failures

### Logging
- **StatusManager** will handle message collection and logging
- Uses `spdlog` package for actual logging
- **Note:** Currently many logging sites have TODO markers during Fornberg implementation

### Smart Pointers
- Use `std::shared_ptr` for shared ownership (domains, methods)
- Use `std::unique_ptr` for exclusive ownership (internal components)
- Prefix member pointers with `mp_` (e.g., `mp_canonical_domain`)

## FFTW & Numerical Methods

### FFTW Coefficient Ordering
- **Forward FFT output order:** `[0, 1, 2, ..., N/2-1, -N/2, -N/2+1, ..., -1]`
- **Normalization:** FFTW does not normalize by default - divide by N for proper DFT coefficients
- **Laurent Series Mapping:** For conformal mapping applications:
  - Positive frequencies (indices 1 to N/2-1) → z^k terms
  - Negative frequencies (indices N/2 to N-1) → z^{-j} terms where j = N-k
  - Index 0 → constant term
  - Index N-1 → z^{-1} coefficient (dominant for conformal maps)

### Fornberg Method Specifics

**Series Representation:**
```
f(z) = Σ a₁,ⱼz^j + Σ_ν Σ_j aν,-j(ρν/(z-cν))^j
```

**Conformal Moduli:**
- c_ν - Centers of circular holes in canonical domain
- ρ_ν - Radii of circular holes in canonical domain
- Total: 3m-6 parameters for m-connected domain (3 degrees of freedom fixed)

**Key Algorithm Components:**
- **Newton iteration** for boundary correspondences S_ν(θ)
- **P matrices** enforce analyticity (zero positive frequency coefficients)
- **CG solver** avoids forming D†D explicitly (memory efficient)
- **Annulus detection** enables specialized m=2 formulation

**Numerical Challenges:**
- System matrices can be ill-conditioned (eigenvalue monitoring)
- Boundary parameters may cluster/reorder during iteration (redistribution)
- Initial guess sensitivity (geometric heuristics help)
- Extreme geometries require careful parameter tuning

### Theodorsen Method Specifics (On Hold)
- **Series form:** f(z) = a₀ + a₋₁z + a₁/z + a₂/z² + ...
- **Numerical stability:** High-order positive powers z^k can cause explosion for |z| > 1
- **Practical evaluation:** Use only constant, linear (z), and first few inverse power (1/z^k) terms
- **Convergence issues:** Extreme ellipse aspect ratios (>2:1) may not converge
- **Current status:** Has formulation error under investigation

## Build & Development Commands

### Building the Project
```bash
# From project root, create and enter build directory
mkdir -p build && cd build

# Configure with cmake (uses clang by default)
CC=clang CXX=clang++ cmake ..

# Build with parallel jobs
make -j9
```

### Running Tests
```bash
# From build directory
./test_exec

# Run specific test suite
./test_exec --gtest_filter="FornbergMCTest.*"
```

### Running GUI
```bash
# From build directory
./conformality_gui
```

### Development Workflow
- All builds happen in the `build/` directory
- Every cmake and make command should be run from the build directory
- Tests are run from the build directory
- Default build type is Release
- Editor: MacVim and shell vim
- CMake automatically generates compile_commands.json for editor integration

### Key Build Targets
- **Main executable:** `conformality_gui` (GUI application)
- **Test executable:** `test_exec` (Google Test suite)

## Documentation & Diagrams

### Design Documentation
- **`design/overview.md`** - Project state and architecture overview
- **`design/roadmap.md`** - Development phases and timeline
- **`design/fornberg/`** - Fornberg method documentation
  - `implementation_plan.md` - Detailed implementation plan (560 lines)
  - `current_status.md` - Current implementation status and next steps

**Important:** The `design/` directory has its own `.git` repository (not a subrepo) for version control during development. This directory and its git history will be removed before publishing the main project.

### Class Diagrams
- Automatic diagrams generated with `clang-uml`
- Configuration: `.clang-uml` in project root
- Output: `docs/diagrams/`
- Claude can generate diagrams defined in config
- **Note:** Diagram generation is currently experimental; clang-uml setup needs refinement

### Code Documentation
- Doxygen-style comments for public APIs
- Implementation notes in source files
- Test files document expected behavior

## Project Management
- **GitHub:** Integration available via `gh` tool
- **Issues:** Can use for planning discussions
- **Milestones:** Track phase completion

## C++ Style Guide

This project uses C++20.

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

## C++ Best Practices

### Memory & Resources
- Use smart pointers (`std::unique_ptr`, `std::shared_ptr`) for dynamic memory
- Prefer `std::vector` over raw arrays (unless compile-time size is known)
- Avoid global variables
- RAII for resource management

### Code Quality
- Always include necessary headers
- Use `const` and `constexpr` extensively (declare const-ness unless value will change)
- Prefer RAII for resource management
- When possible, always go for the best big-O number for performance
- Use Eigen for linear algebra operations
- Use FFTW via FFTWWrapper for FFT operations

### Eigen Usage
- Use `Eigen::MatrixXcd` for complex matrices
- Use `Eigen::VectorXcd` for complex vectors
- Eigen uses value semantics - safe to return from functions
- No need for explicit memory management

## Testing Guidelines

### Framework & Structure
- **Framework:** Google Test (gtest)
- **Location:** All tests in the `tests/` directory
- **Naming:** Test files don't need `test` prefix (since they're in a `tests` directory)
- **Execution:** Run `./test_exec` from the `build/` directory

### Test Coverage Areas

**Currently Tested:**
- **Boundary System:** Component creation, parameterization, sampling
- **Domain Implementations:** Elliptical, circular, and starlike domains
- **Theodorsen Method:** Internal/external mappings (method on hold)
- **Fornberg Components:** Configuration, CGSolver, PMatrixBuilder, FornbergCanonicalDomain
- **Numerical Methods:** FFT wrapper, root finding, coefficient computation
- **Status Management:** Error handling and convergence monitoring

**Needs Additional Testing:**
- **Fornberg Integration:** Full algorithm tests (disabled until implementation complete)
- **Numerical Validation:** Constructive tests validating conformal mapping properties
- **Multiply-Connected Domains:** Various connectivity cases (m=2, 3, 4)

### Test Organization
- Tests are comprehensive and cover all major system components
- Each test file focuses on a specific module or class
- Integration tests verify component interactions
- Numerical tests include tolerance and convergence validation
- Disabled tests marked with `DISABLED_` prefix and TODO comments

### Writing New Tests
- Use descriptive test names: `TEST(ComponentName, TestScenario)`
- Test both success and failure cases
- Use `EXPECT_*` for non-fatal checks, `ASSERT_*` for fatal checks
- Add numerical tolerance checks with `EXPECT_NEAR`
- Document expected behavior in test comments

## Current Development Focus

### Immediate Priorities (Phase 1)

1. **Complete Fornberg Core Methods**
   - Implement `formSystem()`, `solveSystem()`, `newtonUpdate()`
   - Implement boundary sampling and Fourier coefficient computation
   - Complete PMatrixBuilder construction logic

2. **Basic Logging**
   - Add cout-based diagnostics initially
   - StatusManager integration can come in Phase 3

3. **Simple Test Case**
   - Get algorithm running for annulus or simple 2-connected domain
   - Validate that Newton iteration converges

### Next Steps (Phase 2)

1. **Numerical Validation**
   - Create constructive tests that validate algorithm correctness
   - Test convergence properties and stability
   - Validate conformal mapping properties (angle preservation, etc.)
   - Compare algorithm behavior with expected theoretical properties
   - **Note:** We don't have MATLAB license; validation will use constructive tests rather than fixed value comparison

2. **Complete Supporting Methods**
   - Boundary redistribution
   - Inverse mapping
   - Initial guess methods

### Phase 3: Integration & Polish

See `design/roadmap.md` for detailed Phase 3 tasks (logging integration, error handling, code quality).

### Definition of Done (Phase 4)

**Representative multiply-connected domain examples working in GUI**
- Visual verification of correct conformal mapping
- Real-time parameter adjustment
- Example gallery accessible from GUI
- Examples inspired by MATLAB reference repo (but validated constructively)

## Known Issues & Limitations

### Theodorsen Method
- ⚠️ **Critical:** Has formulation error requiring research
- On hold until Fornberg implementation complete

### Fornberg Method
- ~36 TODO/FIXME markers in core implementation
- Needs numerical validation using constructive tests
- Logging integration incomplete
- GUI integration pending

### Platform Support
- Currently macOS only (Metal renderer)
- Cross-platform support planned for future phase

## References

### Academic
- Kropf, E. (2009) "A Fornberg-like Method for the Numerical Conformal Mapping of Bounded Multiply Connected Domains" (Master's Thesis)
- Fornberg, B. (1980) "A Numerical Method for Conformal Mappings"
- Theodorsen, T. (1931) "Theory of Wing Sections of Arbitrary Shape"

### Implementation
- **MATLAB Reference:** `design/fornberg/fornmc/` symlink → `~/Source/gdrive/fornmc/matlab`
- **Design Documentation:** `design/` directory
- **Test Suite:** `tests/` directory

## Contributing & License

**Author:** Everett Kropf (ehkropf@gmail.com)
**License:** AGPL v3 - See LICENSE file
