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

#ifndef FFTW_WRAPPER_HPP
#define FFTW_WRAPPER_HPP

#include <fftw3.h>
#include <vector>
#include <filesystem>

#include "../core/Types.h"

/**
 * @brief A wrapper for the FFTW library that handles 1D FFTs and wisdom management.
 */
class FFTWWrapper
{
private:
    // Wisdom file path
    std::filesystem::path wisdom_file_;

    // FFTW plans for forward and backward transforms
    fftw_plan forward_plan_ = nullptr;
    fftw_plan backward_plan_ = nullptr;

    // Size of the current plan
    size_t current_size_ = 0;

    // Input and output arrays for FFTW
    fftw_complex* input_ = nullptr;
    fftw_complex* output_ = nullptr;

    // Singleton instance
    static FFTWWrapper* instance_;

    /**
     * @brief Load FFTW wisdom from the specified file.
     * @return True if wisdom was successfully loaded, false otherwise.
     */
    bool load_wisdom();

    /**
     * @brief Save FFTW wisdom to the specified file.
     * @return True if wisdom was successfully saved, false otherwise.
     */
    bool save_wisdom();

    /**
     * @brief Create FFTW plans for the given size.
     * @param size Size of the transform.
     */
    void create_plans(size_t size);

    /**
     * @brief Clean up FFTW plans and arrays.
     */
    void cleanup_plans();

    /**
     * @brief Private constructor for singleton pattern.
     * @param wisdom_file Path to the wisdom file.
     */
    FFTWWrapper(const std::filesystem::path& wisdom_file = "fftw.wisdom");

public:
    /**
     * @brief Destructor.
     */
    ~FFTWWrapper();

    /**
     * @brief Get the singleton instance.
     * @param wisdom_file Path to the wisdom file.
     * @return Reference to the singleton instance.
     */
    static FFTWWrapper& get_instance(const std::filesystem::path& wisdom_file = "fftw.wisdom");

    /**
     * @brief Release the singleton instance.
     */
    static void release_instance();

    /**
     * @brief Set the wisdom file path.
     * @param wisdom_file New path to the wisdom file.
     */
    void set_wisdom_file(const std::filesystem::path& wisdom_file);

    /**
     * @brief Get the wisdom file path.
     * @return Current wisdom file path.
     */
    std::filesystem::path get_wisdom_file() const;

    /**
     * @brief Perform a forward FFT (time domain to frequency domain).
     * @param input Input vector of complex values.
     * @return Output vector of complex values.
     */
    std::vector<Complex> forward_fft(const std::vector<Complex>& input);

    /**
     * @brief Perform a backward FFT (frequency domain to time domain).
     * @param input Input vector of complex values.
     * @return Output vector of complex values.
     * @note The output is normalized by dividing by the size of the input vector.
     */
    std::vector<Complex> backward_fft(const std::vector<Complex>& input);

    /**
     * @brief Apply the conjugation operator K to a vector of complex values.
     *
     * The conjugation operator is defined as:
     * K[f(θ)] = -i ∑ σ_m a_m e^(imθ)
     * where a_m are the Fourier coefficients of f and σ_m is the signum function:
     * σ_m = +1 for m > 0, 0 for m = 0, -1 for m < 0.
     *
     * @param input Input vector of complex values.
     * @return Output vector after applying the conjugation operator.
     */
    std::vector<Complex> conjugation_operator(const std::vector<Complex>& input);
};

#endif // FFTW_WRAPPER_HPP
