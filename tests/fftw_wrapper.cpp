/*
 * Copyright (c) 2025, Everett Kropf (ehkropf@gmail.com)
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

#include <gtest/gtest.h>
#include "../src/FFTWWrapper.h"
#include <cmath>
#include <filesystem>

class FFTWWrapperTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        // Use a temporary test-specific wisdom file
        test_wisdom_file_ = std::filesystem::temp_directory_path() / "test_fftw.wisdom";

        // Ensure we start with a fresh instance
        FFTWWrapper::release_instance();
    }

    void TearDown() override
    {
        // Clean up
        FFTWWrapper::release_instance();

        // Remove test wisdom file if it exists
        if (std::filesystem::exists(test_wisdom_file_))
        {
            std::filesystem::remove(test_wisdom_file_);
        }
    }

    std::filesystem::path test_wisdom_file_;
};

TEST_F(FFTWWrapperTest, SingletonBehavior)
{
    auto& instance1 = FFTWWrapper::get_instance(test_wisdom_file_);
    auto& instance2 = FFTWWrapper::get_instance();

    // Both references should point to the same object
    EXPECT_EQ(&instance1, &instance2);

    // The wisdom file path should be set correctly
    EXPECT_EQ(instance1.get_wisdom_file(), test_wisdom_file_);
}

TEST_F(FFTWWrapperTest, WisdomFileHandling)
{
    auto& instance = FFTWWrapper::get_instance(test_wisdom_file_);

    // Test setting a new wisdom file path
    std::filesystem::path new_path = std::filesystem::temp_directory_path() / "new_wisdom.dat";
    instance.set_wisdom_file(new_path);
    EXPECT_EQ(instance.get_wisdom_file(), new_path);

    // Reset to original path for other tests
    instance.set_wisdom_file(test_wisdom_file_);
}

TEST_F(FFTWWrapperTest, ForwardBackwardIdentity)
{
    auto& fftw = FFTWWrapper::get_instance(test_wisdom_file_);

    // Create a simple signal: [1, 0, 0, 0, 0, 0, 0, 0]
    // The FFT should be all ones, and the inverse should return the original
    std::vector<ComplexDouble> signal(8, ComplexDouble(0.0, 0.0));
    signal[0] = ComplexDouble(1.0, 0.0);

    auto forward_result = fftw.forward_fft(signal);

    // All values in the transform should be 1.0
    for (const auto& val : forward_result)
    {
        EXPECT_NEAR(val.real(), 1.0, 1e-10);
        EXPECT_NEAR(val.imag(), 0.0, 1e-10);
    }

    // The backward transform should recover the original signal
    auto backward_result = fftw.backward_fft(forward_result);

    for (size_t i = 0; i < signal.size(); ++i)
    {
        EXPECT_NEAR(backward_result[i].real(), signal[i].real(), 1e-10);
        EXPECT_NEAR(backward_result[i].imag(), signal[i].imag(), 1e-10);
    }
}

TEST_F(FFTWWrapperTest, ForwardTransformSinusoid)
{
    auto& fftw = FFTWWrapper::get_instance(test_wisdom_file_);

    // Create a sinusoidal signal: sin(2π*k/N), k=0,1,...,N-1
    const int N = 16;
    std::vector<ComplexDouble> signal(N);

    for (int k = 0; k < N; ++k)
    {
        double angle = 2.0 * M_PI * k / N;
        signal[k] = ComplexDouble(std::sin(angle), 0.0);
    }

    auto transform = fftw.forward_fft(signal);

    // For a pure sine wave, the FFT should have non-zero components only at
    // frequencies 1 and N-1, and they should be imaginary
    for (int k = 0; k < N; ++k)
    {
        if (k == 1)
        {
            // For sin, the coefficient at frequency 1 should be -0.5i * N
            EXPECT_NEAR(transform[k].real(), 0.0, 1e-10);
            EXPECT_NEAR(transform[k].imag(), -0.5 * N, 1e-10);
        }
        else if (k == N - 1)
        {
            // For sin, the coefficient at frequency N-1 should be 0.5i * N
            EXPECT_NEAR(transform[k].real(), 0.0, 1e-10);
            EXPECT_NEAR(transform[k].imag(), 0.5 * N, 1e-10);
        }
        else
        {
            // All other coefficients should be approximately zero
            EXPECT_NEAR(transform[k].real(), 0.0, 1e-10);
            EXPECT_NEAR(transform[k].imag(), 0.0, 1e-10);
        }
    }
}

TEST_F(FFTWWrapperTest, DifferentSizes)
{
    auto& fftw = FFTWWrapper::get_instance(test_wisdom_file_);

    // Test with size 8
    std::vector<ComplexDouble> signal8(8, ComplexDouble(1.0, 0.0));
    auto result8 = fftw.forward_fft(signal8);
    EXPECT_EQ(result8.size(), 8);

    // Test with size 16
    std::vector<ComplexDouble> signal16(16, ComplexDouble(1.0, 0.0));
    auto result16 = fftw.forward_fft(signal16);
    EXPECT_EQ(result16.size(), 16);

    // Test with size 32
    std::vector<ComplexDouble> signal32(32, ComplexDouble(1.0, 0.0));
    auto result32 = fftw.forward_fft(signal32);
    EXPECT_EQ(result32.size(), 32);
}

TEST_F(FFTWWrapperTest, ConjugationOperator)
{
    auto& fftw = FFTWWrapper::get_instance(test_wisdom_file_);

    // Test the conjugation operator with a simple function
    // f(θ) = cos(θ), whose Hilbert transform is sin(θ)
    const int N = 128;
    std::vector<ComplexDouble> cosine(N);
    std::vector<ComplexDouble> expected_sine(N);

    for (int k = 0; k < N; ++k)
    {
        double angle = 2.0 * M_PI * k / N;
        cosine[k] = ComplexDouble(std::cos(angle), 0.0);
        // The conjugation operator applied to cos(θ) should give sin(θ)
        expected_sine[k] = ComplexDouble(std::sin(angle), 0.0);
    }

    auto result = fftw.conjugation_operator(cosine);

    // Verify the result is close to the expected sine wave
    for (int k = 0; k < N; ++k)
    {
        EXPECT_NEAR(result[k].real(), expected_sine[k].real(), 1e-10);
        EXPECT_NEAR(result[k].imag(), expected_sine[k].imag(), 1e-10);
    }
}

//TEST_F(FFTWWrapperTest, TheodorsenIterationStep)
//{
//    // This test simulates a single step of Theodorsen's method
//    // to verify the conjugation operator works as expected in that context
//
//    auto& fftw = FFTWWrapper::get_instance(test_wisdom_file_);
//
//    const int N = 128;
//    const double a = 1.75;  // Semi-major axis
//    const double b = 1.0;   // Semi-minor axis
//
//    // Initialize arrays
//    std::vector<ComplexDouble> theta(N);
//    std::vector<ComplexDouble> rho(N);
//
//    // Setup initial conditions similar to the Python notebook
//    for (int i = 0; i < N; ++i)
//    {
//        double t = 2.0 * M_PI * i / N;
//        theta[i] = ComplexDouble(t, 0.0);
//
//        // Calculate rho for ellipse
//        double x = a * std::cos(t);
//        double y = b * std::sin(t);
//        double r = std::sqrt(x*x + y*y);
//        rho[i] = ComplexDouble(r, 0.0);
//    }
//
//    // Take log of rho
//    std::vector<ComplexDouble> log_rho(N);
//    for (int i = 0; i < N; ++i)
//    {
//        log_rho[i] = ComplexDouble(std::log(rho[i].real()), 0.0);
//    }
//
//    // Apply conjugation operator (equivalent to K[log(rho)] in Python notebook)
//    auto delta = fftw.conjugation_operator(log_rho);
//
//    // Verify delta[0] is close to 0 (mean of delta should be 0)
//    EXPECT_NEAR(delta[0].real(), 0.0, 1e-10);
//
//    // Verify delta is not all zeros (should have non-zero values for ellipse)
//    bool has_non_zero = false;
//    for (const auto& val : delta)
//    {
//        if (std::abs(val.real()) > 1e-10)
//        {
//            has_non_zero = true;
//            break;
//        }
//    }
//    EXPECT_TRUE(has_non_zero);
//
//    // Calculate phi = delta + theta
//    std::vector<ComplexDouble> phi(N);
//    for (int i = 0; i < N; ++i)
//    {
//        phi[i] = ComplexDouble(delta[i].real() + theta[i].real(), 0.0);
//    }
//
//    // Verify the resulting phi is different from theta
//    for (int i = 0; i < N; ++i)
//    {
//        EXPECT_NE(phi[i].real(), theta[i].real());
//    }
//}
