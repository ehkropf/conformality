# Information for claude

## C++ formatting

These are guidelines, and can be broken if it improves readability.

- Try to keep each line to less than 120 characters.
- Curly braces should always be on their own lines.
- Use spaces around operators `+` and `-` and after commas.
- Use snake_case for variable names.
- Use PascalCase for class names.
- Use camelCase for function names.
- Function definitions should have arguments all on the same line, unless the line would exceed the limit. The closing ')' should always be on the same line as the last argument.
- Continued lines in function definitions should start under the first character of the argment from a previous line. In the case there is no argument after the opening `(` (line space limitations), simply indent the next line by 8 extra spaces instead of 4.
- For the constructor initialisation list, put the colon and the first member init on a separate line to start the list. Each following member initialisation should start the next line with a `,` under the colon followed by the member init.
- Function call arguments should follow the same rules as function definitions.

Example:

```cpp
class MyClass
{
public:
    void myFunction(int param1, int param2)
    {
        // Function body
    }
};
```

## C++ best practices

- Always include necessary headers.
- Use `const` and `constexpr` where appropriate.
- Prefer `std::vector` over raw arrays.
- Use smart pointers (`std::unique_ptr`, `std::shared_ptr`) for dynamic memory management.
- Avoid global variables.
