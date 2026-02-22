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

#include <stdexcept>

namespace conformality::examples
{

ThesisExamplePreset ThesisExamples::getExample(int /*exampleNumber*/)
{
    throw std::runtime_error("ThesisExamples::getExample not implemented");
}

std::vector<int> ThesisExamples::availableExamples()
{
    throw std::runtime_error("ThesisExamples::availableExamples not implemented");
}

ThesisExamplePreset ThesisExamples::makeExample2()
{
    throw std::runtime_error("ThesisExamples::makeExample2 not implemented");
}

ThesisExamplePreset ThesisExamples::makeExample3()
{
    throw std::runtime_error("ThesisExamples::makeExample3 not implemented");
}

ThesisExamplePreset ThesisExamples::makeExample4()
{
    throw std::runtime_error("ThesisExamples::makeExample4 not implemented");
}

ThesisExamplePreset ThesisExamples::makeExample5()
{
    throw std::runtime_error("ThesisExamples::makeExample5 not implemented");
}

FornbergMCConfiguration ThesisExamples::makeConfig(int /*N*/)
{
    throw std::runtime_error("ThesisExamples::makeConfig not implemented");
}

} // namespace conformality::examples
