# Thesis Example Preset System Design

**Issue:** #81
**Date:** 2026-02-22

## Summary

Create a preset system providing hardcoded thesis example domain configurations, shared between GUI and CLI. Includes a new `InvertedEllipseComponent` boundary type and shared boundary factory helpers.

## File Structure

```
src/domains/
  InvertedEllipseComponent.h     # New boundary type (binvellip)
  InvertedEllipseComponent.cpp

src/examples/
  ThesisExamples.h               # Preset struct + static class
  ThesisExamples.cpp
  BoundaryHelpers.h              # Shared boundary factory functions
  BoundaryHelpers.cpp
```

## InvertedEllipseComponent

New `BoundaryComponent` subclass in `src/domains/`. Implements the MATLAB `binvellip` parameterization.

```cpp
class InvertedEllipseComponent : public BoundaryComponent
{
public:
    InvertedEllipseComponent(Complex center, double alpha, double rotation = 0.0);

    Complex evaluate(double t) const override;
    Complex evaluateDerivative(double t) const override;
    std::vector<Complex> sample(size_t numPoints) const override;
    double findParameterization(const Complex& z) const override;

private:
    Complex m_center;
    double m_alpha;
    double m_rotation;
};
```

**Math (rotation = 0):**
```
z(t)  = 1 / (alpha * cos(t) - i * sin(t)) + center
z'(t) = (alpha * sin(t) + i * cos(t)) / (alpha * cos(t) - i * sin(t))^2
```

**General rotation R:**
```
denom = alpha*cos(R)*cos(t) - sin(R)*sin(t) - i*(cos(R)*sin(t) + alpha*sin(R)*cos(t))
z(t)  = 1/denom + center
z'(t) = (alpha*cos(R)*sin(t) + sin(R)*cos(t) + i*(cos(R)*cos(t) - alpha*sin(R)*sin(t))) / denom^2
```

`sample()` uses uniform spacing. `findParameterization()` uses `RootFinder::ternarySearch`.

## BoundaryHelpers

```cpp
namespace conformality::examples
{

std::shared_ptr<Boundary> createCircularBoundary(Complex center, double radius);

std::shared_ptr<Boundary> createEllipseBoundary(
    Complex center, double semi_major, double semi_minor, double rotation = 0.0);

std::shared_ptr<Boundary> createInvertedEllipseBoundary(
    Complex center, double alpha, double rotation = 0.0);

} // namespace conformality::examples
```

- `createCircularBoundary` and `createEllipseBoundary` wrap `AnalyticBoundaryComponent` with lambdas
- `createInvertedEllipseBoundary` wraps `InvertedEllipseComponent`
- All return `shared_ptr<Boundary>` with a single component
- Eliminates duplication of `createCircularBoundary` across test files

## ThesisExamples

```cpp
namespace conformality::examples
{

struct ThesisExamplePreset
{
    std::string name;
    std::string description;
    std::shared_ptr<MultiplyConnectedDomain> target_domain;
    std::vector<Complex> initial_centers;
    std::vector<double> initial_radii;
    FornbergMCConfiguration config;
};

class ThesisExamples
{
public:
    static ThesisExamplePreset getExample(int exampleNumber);
    static std::vector<int> availableExamples();  // {3, 5, 2, 4}
};

} // namespace conformality::examples
```

### Presets

| n | Name | m | Outer | Inner | N |
|---|------|---|-------|-------|---|
| 3 | Identity (m=4) | 4 | Circle(0, 1) | 3 circles | 256 |
| 5 | Ellipses (m=3) | 3 | Ellipse(0, 2, 1.5) | 2 ellipses | 256 |
| 2 | Mixed (m=4) | 4 | InvEllipse(0, 0.3) | 3 ellipses | 128 |
| 4 | High connectivity (m=7) | 7 | Ellipse(0, 2, 1) | 6 ellipses | 128 |

`availableExamples()` returns `{3, 5, 2, 4}` — ordered simplest-first for GUI menu presentation.

### Preset Parameters (from MATLAB th_gen_ex*.m)

**Ex 3 (Identity, m=4):** All circles. Outer: Circle(0, 1). Inner: Circle(-0.5, 0.25), Circle(0.25+0.43i, 0.25), Circle(0.25-0.43i, 0.25). Initial guesses offset ~0.1 from true centers: (-0.4, 0.25), (0.35+0.43i, 0.25), (0.35-0.43i, 0.25).

**Ex 5 (Ellipses, m=3):** Outer: Ellipse(0, 2, 1.5). Inner: Ellipse(-0.8, 3/16, 3/8, 0), Ellipse(0.7-0.1i, 3/8, 3/16, pi/4). Initial guesses: (-0.47, 0.15), (0.4-0.1i, 0.15).

**Ex 2 (Mixed, m=4):** Outer: InvEllipse(0, 0.3). Inner: Ellipse(1+0.3i, 3/4, 3/8, pi/4), Ellipse(1.7-0.7i, 1/2, 1/4, pi/4), Ellipse(-1.7, 3/8, 3/4, 0). Initial guesses: (0.6+0.1i, 0.14), (0.77-0.2i, 0.05), (-0.7, 0.2).

**Ex 4 (High connectivity, m=7):** Outer: Ellipse(0, 2, 1). Six inner ellipses with various rotations. See MATLAB `th_gen_ex4.m` for full parameters.

### Configuration (all presets)

- `newton_tolerance = 1e-14`
- `cgm_tolerance = 1e-15`
- `normalization_condition = {1, 0, 0}`
- Other values use `FornbergMCConfiguration` defaults

### Error Handling

`getExample()` throws `std::invalid_argument` for unsupported example numbers.

## CMake Integration

- `InvertedEllipseComponent` added to `domain_objects`
- New `example_objects` object library in `src/examples/`
- `example_objects` links against `domain_objects`, `method_objects`
- GUI target links against `example_objects`

## Testing

- **InvertedEllipseComponent:** evaluate at known angles, derivative vs finite difference, sample count, findParameterization round-trip
- **BoundaryHelpers:** circle radius correctness, ellipse semi-axes at cardinal points
- **ThesisExamples:** each preset constructs successfully, correct connectivity, initial guess count matches inner boundary count, config values match expected

## Not In Scope

- Running `FornbergMC::compute()` on presets
- GUI modifications (#82, #83)
- CLI program (#84)
- Updating existing tests to use BoundaryHelpers (cleanup, not required)
