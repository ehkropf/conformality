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

#pragma once

#include <iomanip>
#include <sstream>
#include <string>
#include <variant>
#include <vector>

using MethodInfoValue = std::variant<int, double, bool, std::string>;

struct MethodInfoField
{
    std::string label;
    MethodInfoValue value;
};

struct MethodInfo
{
    std::string name;
    std::vector<MethodInfoField> parameters;
    std::vector<MethodInfoField> results;
};

template <class... Ts>
struct overloaded : Ts... { using Ts::operator()...; };

inline std::string formatMethodInfoValue(const MethodInfoValue& value)
{
    return std::visit(overloaded{
        [](int v) { return std::to_string(v); },
        [](double v) { std::ostringstream s; s << std::scientific << std::setprecision(2) << v; return s.str(); },
        [](bool v) -> std::string { return v ? "Yes" : "No"; },
        [](const std::string& v) { return v; }
    }, value);
}
