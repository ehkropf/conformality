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

#include "FFTWWrapper.h"

// Initialize static member
FFTWWrapper* FFTWWrapper::instance_ = nullptr;

// Private constructor implementation
FFTWWrapper::FFTWWrapper(const std::filesystem::path& wisdom_file)
    : wisdom_file_(wisdom_file)
{
    // Initialize FFTW threads (optional)
    // fftw_init_threads();
    // fftw_plan_with_nthreads(omp_get_max_threads());
}

// Destructor implementation
FFTWWrapper::~FFTWWrapper()
{
    cleanup_plans();
    // Cleanup FFTW threads (if used)
    // fftw_cleanup_threads();
    fftw_cleanup();
}

// Singleton instance getter
FFTWWrapper& FFTWWrapper::get_instance(const std::filesystem::path& wisdom_file)
{
    if (!instance_)
    {
        instance_ = new FFTWWrapper(wisdom_file);
    }
    return *instance_;
}

// Release singleton instance
void FFTWWrapper::release_instance()
{
    if (instance_)
    {
        delete instance_;
        instance_ = nullptr;
    }
}

// Wisdom file getter and setter
void FFTWWrapper::set_wisdom_file(const std::filesystem::path& wisdom_file)
{
    wisdom_file_ = wisdom_file;
}

std::filesystem::path FFTWWrapper::get_wisdom_file() const
{
    return wisdom_file_;
}

// Wisdom management methods
bool FFTWWrapper::load_wisdom()
{
    if (std::filesystem::exists(wisdom_file_))
    {
        int result = fftw_import_wisdom_from_filename(wisdom_file_.c_str());
        return result != 0;
    }
    return false;
}

bool FFTWWrapper::save_wisdom()
{
    int result = fftw_export_wisdom_to_filename(wisdom_file_.c_str());
    return result != 0;
}

// Plan management methods
void FFTWWrapper::create_plans(size_t size)
{
    // Clean up old plans if they exist
    cleanup_plans();

    // Allocate memory for input and output arrays
    input_ = fftw_alloc_complex(size);
    output_ = fftw_alloc_complex(size);

    if (!input_ || !output_)
    {
        throw std::runtime_error("Failed to allocate memory for FFTW arrays");
    }

    // Load existing wisdom
    load_wisdom();

    // Create new plans using FFTW_MEASURE for better performance
    forward_plan_ = fftw_plan_dft_1d(size, input_, output_, FFTW_FORWARD, FFTW_MEASURE);
    backward_plan_ = fftw_plan_dft_1d(size, output_, input_, FFTW_BACKWARD, FFTW_MEASURE);

    if (!forward_plan_ || !backward_plan_)
    {
        cleanup_plans();
        throw std::runtime_error("Failed to create FFTW plans");
    }

    // Save the new wisdom
    save_wisdom();

    // Update current size
    current_size_ = size;
}

void FFTWWrapper::cleanup_plans()
{
    if (forward_plan_)
    {
        fftw_destroy_plan(forward_plan_);
        forward_plan_ = nullptr;
    }

    if (backward_plan_)
    {
        fftw_destroy_plan(backward_plan_);
        backward_plan_ = nullptr;
    }

    if (input_)
    {
        fftw_free(input_);
        input_ = nullptr;
    }

    if (output_)
    {
        fftw_free(output_);
        output_ = nullptr;
    }

    current_size_ = 0;
}

// FFT methods
std::vector<ComplexDouble> FFTWWrapper::forward_fft(const std::vector<ComplexDouble>& input)
{
    size_t size = input.size();

    // Ensure we have plans for this size
    if (size != current_size_)
    {
        create_plans(size);
    }

    // Copy input data to FFTW input array
    for (size_t i = 0; i < size; ++i)
    {
        input_[i][0] = input[i].real();
        input_[i][1] = input[i].imag();
    }

    // Execute the forward plan
    fftw_execute(forward_plan_);

    // Copy output data to result vector
    std::vector<ComplexDouble> result(size);
    for (size_t i = 0; i < size; ++i)
    {
        result[i] = ComplexDouble(output_[i][0], output_[i][1]);
    }

    return result;
}

std::vector<ComplexDouble> FFTWWrapper::backward_fft(const std::vector<ComplexDouble>& input)
{
    size_t size = input.size();

    // Ensure we have plans for this size
    if (size != current_size_)
    {
        create_plans(size);
    }

    // Copy input data to FFTW output array (since backward transform starts from output)
    for (size_t i = 0; i < size; ++i)
    {
        output_[i][0] = input[i].real();
        output_[i][1] = input[i].imag();
    }

    // Execute the backward plan
    fftw_execute(backward_plan_);

    // Copy input data (which now contains the result) to result vector
    std::vector<ComplexDouble> result(size);
    for (size_t i = 0; i < size; ++i)
    {
        // Normalize by dividing by size
        result[i] = ComplexDouble(input_[i][0] / size, input_[i][1] / size);
    }

    return result;
}

// Conjugation operator
std::vector<ComplexDouble> FFTWWrapper::conjugation_operator(const std::vector<ComplexDouble>& input)
{
    // First compute the FFT
    std::vector<ComplexDouble> fourier_coeffs = forward_fft(input);
    size_t size = fourier_coeffs.size();
    size_t half_size = size / 2;

    // Apply the signum function to the Fourier coefficients
    // Set the DC component (m = 0) to zero
    fourier_coeffs[0] = ComplexDouble(0.0, 0.0);

    // For m > 0, multiply by -i
    for (size_t i = 1; i < half_size; ++i)
    {
        fourier_coeffs[i] = fourier_coeffs[i] * ComplexDouble(0.0, -1.0);
    }

    // For m < 0, multiply by i
    for (size_t i = half_size; i < size; ++i)
    {
        fourier_coeffs[i] = fourier_coeffs[i] * ComplexDouble(0.0, 1.0);
    }

    // Compute the inverse FFT
    return backward_fft(fourier_coeffs);
}
