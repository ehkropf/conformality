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

#include "../examples/ThesisExamples.h"
#include "../domains/FornbergCanonicalDomain.h"
#include "../methods/FornbergMC.h"
#include "../methods/PMatrixBuilder.h"
#include "../numerics/CGSolver.h"
#include "../core/ConformalMap.h"
#include "../core/StatusManager.h"

#include <iostream>
#include <map>
#include <stdexcept>
#include <string>

using namespace conformality::examples;

namespace
{

const std::map<std::string, int> EXAMPLE_NAMES = {
    {"thesis1", 1},
    {"thesis2", 2},
    {"thesis3", 3},
    {"thesis4", 4},
    {"thesis5", 5},
};

void printUsage(const char* program_name)
{
    std::cout << "Usage: " << program_name << " [OPTIONS]\n"
              << "\n"
              << "Conformality CLI - run FornbergMC conformal mapping on thesis examples.\n"
              << "\n"
              << "Options:\n"
              << "  --example <name>   Run a specific example (default: thesis3)\n"
              << "  --list-examples    List available examples\n"
              << "  --help             Show this help message\n"
              << "\n"
              << "Examples:\n"
              << "  " << program_name << "                        # Run default (thesis3)\n"
              << "  " << program_name << " --example thesis5      # Run thesis example 5\n"
              << "  " << program_name << " --list-examples        # Show available examples\n";
}

void printAvailableExamples()
{
    std::cout << "Available examples:\n";
    for (int num : ThesisExamples::availableExamples())
    {
        auto preset = ThesisExamples::getExample(num);
        std::cout << "  thesis" << num << "  - " << preset.description << "\n";
    }
}

int runExample(const std::string& example_name)
{
    auto it = EXAMPLE_NAMES.find(example_name);
    if (it == EXAMPLE_NAMES.end())
    {
        std::cerr << "Error: unknown example '" << example_name << "'\n";
        std::cerr << "Use --list-examples to see available examples.\n";
        return 1;
    }

    int example_number = it->second;
    auto preset = ThesisExamples::getExample(example_number);

    std::cout << "=== " << preset.name << " ===\n";
    std::cout << "Description: " << preset.description << "\n";
    std::cout << "Connectivity: " << preset.target_domain->getConnectivity() << "\n";
    std::cout << "N (points per boundary): " << preset.config.N << "\n";
    std::cout << "\n";

    // Create canonical domain from initial guesses
    auto source_domain = std::make_shared<FornbergCanonicalDomain>(
        preset.initial_centers, preset.initial_radii, preset.config.N);

    // Create method with status manager for console logging
    auto method = std::make_shared<FornbergMC>(preset.config);
    auto status_manager = std::make_shared<StatusManager>();
    status_manager->enableLogging(LogOutput::CONSOLE);
    method->setStatusManager(status_manager);

    // Create and compute the conformal map
    ConformalMap map(source_domain, preset.target_domain, method);

    std::cout << "Computing conformal map...\n\n";
    map.compute();

    // Print results
    std::cout << "\n=== Results ===\n";
    std::cout << "Converged: " << (method->hasConverged() ? "yes" : "no") << "\n";
    std::cout << "Iterations: " << method->getResidualHistory().size() << "\n";
    std::cout << "Final residual: " << method->getCurrentResidual() << "\n";

    // Print conformal moduli
    const auto& moduli = method->getConformalModuli();
    std::cout << "\nConformal moduli (interleaved [c, rho] for each inner boundary):\n";
    for (Eigen::Index i = 0; i < moduli.size(); i += 2)
    {
        int boundary_idx = static_cast<int>(i / 2) + 2; // inner boundaries start at index 2
        Complex center = moduli(i);
        double radius = moduli(i + 1).real();
        std::cout << "  Boundary " << boundary_idx << ": c = ("
                  << center.real() << ", " << center.imag() << "), rho = " << radius << "\n";
    }

    return method->hasConverged() ? 0 : 1;
}

} // anonymous namespace

int main(int argc, char* argv[])
{
    std::string example_name = "thesis3";

    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];

        if (arg == "--help" || arg == "-h")
        {
            printUsage(argv[0]);
            return 0;
        }
        else if (arg == "--list-examples")
        {
            printAvailableExamples();
            return 0;
        }
        else if (arg == "--example")
        {
            if (i + 1 >= argc)
            {
                std::cerr << "Error: --example requires a value\n";
                return 1;
            }
            example_name = argv[++i];
        }
        else
        {
            std::cerr << "Error: unknown option '" << arg << "'\n";
            std::cerr << "Use --help for usage information.\n";
            return 1;
        }
    }

    try
    {
        return runExample(example_name);
    }
    catch (const std::invalid_argument& e)
    {
        std::cerr << "Validation error: " << e.what() << "\n";
        return 2;
    }
    catch (const std::runtime_error& e)
    {
        std::cerr << "Runtime error: " << e.what() << "\n";
        return 1;
    }
}
