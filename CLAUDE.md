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
- Test files don't need `test` in the name, since there in a `tests` directory
