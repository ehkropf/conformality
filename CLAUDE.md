# Claude Project Guide

## Project Status & Documentation

**Current Phase:** Phase 1
**Project Overview:** `design/high_level_overview.md`
**Phase 1 Documentation:** `design/phase1/` directory
**Map Design Discussion:** `design/phase1/map_hierarchy_design_discussion.md` (led to major refactor)

## Project architecture

TBD coherently, but see design directory.

## Dependencies & Libraries

- **Third-party included libraries**
  - imgui
  - implot
  - spdlog
- **Libraries that must be pre-installed**
  - FFTW (needs compiled for target system)

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

## Build & development commands
- Builds will be done in a `build` directory
  - Every cmake and make command should be given from this directory
  - Tests will be run from this directory
- Use `cmake` to configure the build
  - CC=clang
  - CXX=clang++
  - Release is default build type
- Use `make -j9` to build
- No IDE used, just MacVim and shell vim.

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

## C++ Best Practices

### Memory & Resources

- Use smart pointers (`std::unique_ptr`, `std::shared_ptr`) for dynamic memory
- Prefer `std::vector` over raw arrays (unless compile-time size is known)
- Avoid global variables

### Code Quality

- Always include necessary headers
- Use `const` and `constexpr` extensively (declare const-ness unless value will change)
- Prefer RAII for resource management
- When possible, always go for the best big-O number for performance

## Testing guidelines

- Using gtest
- All tests go in the `tests` directory.
- Test files don't need `test` in the name, since they are in a `tests` directory
