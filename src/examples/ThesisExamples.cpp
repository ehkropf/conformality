/*
 * Copyright © 2025, Everett Kropf (ehkropf@gmail.com)
 *
 * This file is part of Conformality.
 * Conformality is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Conformality is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Affero General Public License along
 * with Conformality. If not, see <https://www.gnu.org/licenses/>.
 */

#include "ThesisExamples.h"
#include "BoundaryHelpers.h"

#include <stdexcept>

namespace conformality::examples
{

ThesisExamplePreset ThesisExamples::getExample(int exampleNumber)
{
    switch (exampleNumber)
    {
        case 1: return makeExample1();
        case 2: return makeExample2();
        case 3: return makeExample3();
        case 4: return makeExample4();
        case 5: return makeExample5();
        default:
            throw std::invalid_argument(
                "ThesisExamples: unsupported example number " + std::to_string(exampleNumber)
                + " (available: 1, 2, 3, 4, 5)");
    }
}

std::vector<int> ThesisExamples::availableExamples()
{
    return {1, 3, 5, 2, 4};
}

// MATLAB th_gen_ex1.m: m=3, spline outer boundary + two inner ellipses
ThesisExamplePreset ThesisExamples::makeExample1()
{
    // Outer boundary: periodic cubic spline through 11 control points
    // Control points extracted from th_gen_ex1_spline.mat
    std::vector<double> xpts = {
        1.956140, 1.570175, 0.710526, 0.008772, -0.412281, -1.289474,
        -1.798246, -2.026316, -1.149123, 0.692982, 1.728070, 1.956140
    };
    std::vector<double> ypts = {
        0.043860, 0.500000, 0.657895, 0.815789, 1.429825, 1.605263,
        0.710526, -0.622807, -1.675439, -1.763158, -1.061404, 0.043860
    };
    auto outer = createSplineBoundary(xpts, ypts, 256);

    // C2 = bellipse([-.9  3/16   3/4   0], N, [-.57-.1i, .2])
    auto inner1 = createEllipseBoundary(Complex(-0.9, 0), 3.0 / 16.0, 3.0 / 4.0, 0.0);

    // C3 = bellipse([.7-.5i   3/4   3/16  pi/4], N, [.25-.5i, .22])
    auto inner2 = createEllipseBoundary(Complex(0.7, -0.5), 3.0 / 4.0, 3.0 / 16.0, M_PI / 4.0);

    auto domain = std::make_shared<MultiplyConnectedDomain>(
        std::vector<std::shared_ptr<Boundary>>{outer, inner1, inner2});

    return ThesisExamplePreset{
        "Thesis Example 1",
        "Spline outer boundary (m=3): periodic cubic spline with two inner ellipses",
        domain,
        {Complex(-0.57, -0.1), Complex(0.25, -0.5)},
        {0.2, 0.22},
        makeConfig(256)
    };
}

// MATLAB th_gen_ex2.m: m=4, mixed (inverted ellipse outer + inner ellipses)
ThesisExamplePreset ThesisExamples::makeExample2()
{
    // C1 = binvellip([0  .3], N)
    auto outer = createInvertedEllipseBoundary(Complex(0, 0), 0.3);

    // C2 = bellipse([1+.3i   3/4   3/8  pi/4], N, [.6+.1i, .14])
    auto inner1 = createEllipseBoundary(Complex(1.0, 0.3), 3.0 / 4.0, 3.0 / 8.0, M_PI / 4.0);

    // C3 = bellipse([1.7-.7i   1/2   1/4  pi/4], N, [.77-.2i, .05])
    auto inner2 = createEllipseBoundary(Complex(1.7, -0.7), 1.0 / 2.0, 1.0 / 4.0, M_PI / 4.0);

    // C4 = bellipse([-1.7   3/8   3/4], N, [-.7, .2])
    auto inner3 = createEllipseBoundary(Complex(-1.7, 0), 3.0 / 8.0, 3.0 / 4.0);

    auto domain = std::make_shared<MultiplyConnectedDomain>(
        std::vector<std::shared_ptr<Boundary>>{outer, inner1, inner2, inner3});

    return ThesisExamplePreset{
        "Thesis Example 2",
        "Mixed boundaries (m=4): inverted ellipse outer with three inner ellipses",
        domain,
        {Complex(0.6, 0.1), Complex(0.77, -0.2), Complex(-0.7, 0)},
        {0.14, 0.05, 0.2},
        makeConfig(256)
    };
}

// MATLAB th_gen_ex3.m: m=4, identity map (all circles)
ThesisExamplePreset ThesisExamples::makeExample3()
{
    // C1 = bcircle([0  1], N)
    auto outer = createCircularBoundary(Complex(0, 0), 1.0);

    // C2 = bcircle([-.5   .25], N, [-.4  .25])
    auto inner1 = createCircularBoundary(Complex(-0.5, 0), 0.25);

    // C3 = bcircle([.25+.43i  .25], N, [.35+.43i  .25])
    auto inner2 = createCircularBoundary(Complex(0.25, 0.43), 0.25);

    // C4 = bcircle([.25-.43i  .25], N, [.35-.43i  .25])
    auto inner3 = createCircularBoundary(Complex(0.25, -0.43), 0.25);

    auto domain = std::make_shared<MultiplyConnectedDomain>(
        std::vector<std::shared_ptr<Boundary>>{outer, inner1, inner2, inner3});

    return ThesisExamplePreset{
        "Thesis Example 3",
        "Identity map (m=4): unit circle outer with three inner circles",
        domain,
        {Complex(-0.4, 0), Complex(0.35, 0.43), Complex(0.35, -0.43)},
        {0.25, 0.25, 0.25},
        makeConfig(256)
    };
}

// MATLAB th_gen_ex4.m: m=7, high connectivity (ellipses)
// Note: C8 is defined in MATLAB but not used in tregion (C1-C7 only)
ThesisExamplePreset ThesisExamples::makeExample4()
{
    // C1 = bellipse([0  2   1], N)
    auto outer = createEllipseBoundary(Complex(0, 0), 2.0, 1.0);

    // C2 = bellipse([1.2+.3i  1/4  1/8  -pi/12], N, [.8+.15i, .1])
    auto inner1 = createEllipseBoundary(Complex(1.2, 0.3), 1.0 / 4.0, 1.0 / 8.0, -M_PI / 12.0);

    // C3 = bellipse([1-.3i  1/4  1/8  pi/12], N, [.7-.15i, .1])
    auto inner2 = createEllipseBoundary(Complex(1.0, -0.3), 1.0 / 4.0, 1.0 / 8.0, M_PI / 12.0);

    // C4 = bellipse([.5   1/4  1/8], N, [.4,  .1])
    auto inner3 = createEllipseBoundary(Complex(0.5, 0), 1.0 / 4.0, 1.0 / 8.0);

    // C5 = bellipse([-.8  1/8  1/4], N, [-.6,  .1])
    auto inner4 = createEllipseBoundary(Complex(-0.8, 0), 1.0 / 8.0, 1.0 / 4.0);

    // C6 = bellipse([-.25+.5i  1/4  1/8], N, [-.2+.4i,  .15])
    auto inner5 = createEllipseBoundary(Complex(-0.25, 0.5), 1.0 / 4.0, 1.0 / 8.0);

    // C7 = bellipse([-.3-.5i  1/4  1/8  pi/6], N, [-.3-.4i,  .15])
    auto inner6 = createEllipseBoundary(Complex(-0.3, -0.5), 1.0 / 4.0, 1.0 / 8.0, M_PI / 6.0);

    auto domain = std::make_shared<MultiplyConnectedDomain>(
        std::vector<std::shared_ptr<Boundary>>{outer, inner1, inner2, inner3, inner4, inner5, inner6});

    return ThesisExamplePreset{
        "Thesis Example 4",
        "High connectivity (m=7): ellipse outer with six inner ellipses",
        domain,
        {Complex(0.8, 0.15), Complex(0.7, -0.15), Complex(0.4, 0),
         Complex(-0.6, 0), Complex(-0.2, 0.4), Complex(-0.3, -0.4)},
        {0.1, 0.1, 0.1, 0.1, 0.15, 0.15},
        makeConfig(128)
    };
}

// MATLAB th_gen_ex5.m: m=3, simple ellipses
ThesisExamplePreset ThesisExamples::makeExample5()
{
    // C1 = bellipse([0  2   3/2], N)
    auto outer = createEllipseBoundary(Complex(0, 0), 2.0, 3.0 / 2.0);

    // C2 = bellipse([-.8  3/16   3/8   0], N, [-.47  .15])
    auto inner1 = createEllipseBoundary(Complex(-0.8, 0), 3.0 / 16.0, 3.0 / 8.0, 0.0);

    // C3 = bellipse([.7-.1i   3/8   3/16  pi/4], N, [.4-.1i  .15])
    auto inner2 = createEllipseBoundary(Complex(0.7, -0.1), 3.0 / 8.0, 3.0 / 16.0, M_PI / 4.0);

    auto domain = std::make_shared<MultiplyConnectedDomain>(
        std::vector<std::shared_ptr<Boundary>>{outer, inner1, inner2});

    return ThesisExamplePreset{
        "Thesis Example 5",
        "Ellipses (m=3): ellipse outer with two inner ellipses",
        domain,
        {Complex(-0.47, 0), Complex(0.4, -0.1)},
        {0.15, 0.15},
        makeConfig(256)
    };
}

FornbergMCConfiguration ThesisExamples::makeConfig(int N)
{
    FornbergMCConfiguration config;
    config.N = N;
    config.newton_tolerance = 1e-14;
    config.cgm_tolerance = 1e-15;
    config.initial_guess_method = FornbergMCConfiguration::InitialGuessMethod::MANUAL;
    return config;
}

} // namespace conformality::examples
