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

#include <Eigen/Dense>
#include <nlohmann/json.hpp>
#include <complex>
#include <fstream>
#include <stdexcept>
#include <string>

/**
 * @brief Loads MATLAB/Octave reference data from JSON files for comparison tests.
 *
 * JSON format conventions:
 * - Complex numbers: [real, imag]
 * - Matrices: array of row arrays
 * - Vectors: array of [real, imag] pairs
 */
class ReferenceDataLoader
{
public:
    explicit ReferenceDataLoader(const std::string& json_path)
    {
        std::ifstream file(json_path);
        if (!file.is_open())
        {
            throw std::runtime_error("Cannot open reference data file: " + json_path);
        }
        try
        {
            m_data = nlohmann::json::parse(file);
        }
        catch (const nlohmann::json::parse_error& e)
        {
            throw std::runtime_error("Failed to parse reference data file: " + json_path + " - " + e.what());
        }
    }

    Eigen::MatrixXcd getComplexMatrix(const std::string& key) const
    {
        const auto& mat = m_data.at("matrices").at(key);
        int rows = static_cast<int>(mat.size());
        if (rows == 0)
        {
            return Eigen::MatrixXcd(0, 0);
        }
        int cols = static_cast<int>(mat[0].size());
        Eigen::MatrixXcd result(rows, cols);

        for (int r = 0; r < rows; ++r)
        {
            for (int c = 0; c < cols; ++c)
            {
                const auto& val = mat[r][c];
                result(r, c) = std::complex<double>(val[0].get<double>(), val[1].get<double>());
            }
        }
        return result;
    }

    Eigen::MatrixXd getRealMatrixFromComplex(const std::string& key) const
    {
        auto cmat = getComplexMatrix(key);
        return cmat.real();
    }

    Eigen::VectorXcd getComplexVector(const std::string& key) const
    {
        const auto& vec = m_data.at("vectors").at(key);
        int n = static_cast<int>(vec.size());
        Eigen::VectorXcd result(n);
        for (int i = 0; i < n; ++i)
        {
            const auto& val = vec[i];
            result(i) = std::complex<double>(val[0].get<double>(), val[1].get<double>());
        }
        return result;
    }

    Eigen::VectorXd getRealVectorFromComplex(const std::string& key) const
    {
        auto cvec = getComplexVector(key);
        return cvec.real();
    }

    int getMetadataInt(const std::string& key) const
    {
        return m_data.at("metadata").at(key).get<int>();
    }

    double getMetadataDouble(const std::string& key) const
    {
        return m_data.at("metadata").at(key).get<double>();
    }

    std::string getMetadataString(const std::string& key) const
    {
        return m_data.at("metadata").at(key).get<std::string>();
    }

    std::vector<double> getMetadataArray(const std::string& key) const
    {
        return m_data.at("metadata").at(key).get<std::vector<double>>();
    }

    double getScalar(const std::string& key) const
    {
        return m_data.at("scalars").at(key).get<double>();
    }

    bool hasSection(const std::string& section) const
    {
        return m_data.contains(section);
    }

    bool hasKey(const std::string& section, const std::string& key) const
    {
        return m_data.contains(section) && m_data.at(section).contains(key);
    }

private:
    nlohmann::json m_data;
};
