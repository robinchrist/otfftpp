/******************************************************************************
*  OTFFT Header Version 11.4xv
*
*  Copyright (c) 2015 OK Ojisan(Takuya OKAHISA)
*  Released under the MIT license
*  http://opensource.org/licenses/mit-license.php
******************************************************************************/

#ifndef otfft_h
#define otfft_h

#include "otfftpp/otfft_fwd.h"
#include "otfftpp/otfft_complex.h"

#include <memory>
#include <cstdlib>
#include <cstddef>
#include <cstdint>

namespace OTFFT
{
    enum class TransformationType
    {
        // Complex Discrete Fourier Transform
        TRANSFORM_FFT_COMPLEX,
        // Real Discrete Fourier Transform
        TRANSFORM_FFT_REAL,
        // Discrete Cosine Transform (DCT-II)
        TRANSFORM_DCT,
        // Bluestein's FFT
        TRANSFORM_BLUESTEIN
    };

    /**
     * @brief Complex In-Place Discrete Fourier Transform
     */
    class ComplexFFT
    {
    public:
        virtual ~ComplexFFT() = default;

        /**
         * @brief Setup the sequence length of the algorithm (up to 2^30)
         */
        virtual void setup(int n) = 0;

        /**
         * @brief Discrete forward transformation (with normalization by 1/N)
         */
        virtual void fwd(complex_vector  x) const noexcept = 0;
        /**
         * @brief Discrete forward transformation (without normalization)
         */
        virtual void fwd0(complex_vector x) const noexcept = 0;
        /**
         * @brief Discrete forward transformation (unitary transformation, with normalization by sqrt(1/N))
         */
        virtual void fwdu(complex_vector x) const noexcept = 0;
        /**
         * @brief Discrete forward transformation (with normalization by 1/N)
         */
        virtual void fwdn(complex_vector x) const noexcept = 0;
        /**
         * @brief Discrete inverse transformation (without normalization)
         */
        virtual void inv(complex_vector  x) const noexcept = 0;
        /**
         * @brief Discrete inverse transformation (without normalization)
         */
        virtual void inv0(complex_vector x) const noexcept = 0;
        /**
         * @brief Discrete inverse transformation (unitary transformation)
         */
        virtual void invu(complex_vector x) const noexcept = 0;
        /**
         * @brief Discrete inverse transformation (with normalization by 1/N)
         */
        virtual void invn(complex_vector x) const noexcept = 0;
    };

    /**
     * @brief Real Out-Of-Place Discrete Fourier Transform
     */
    class RealFFT
    {
    public:
        virtual ~RealFFT() = default;

        /**
         * @brief Setup the sequence length of the algorithm (up to 2^30)
         */
        virtual void setup(int n) = 0;

        /**
         * @brief Discrete forward transformation (with normalization by 1/N)
         */
        virtual void fwd(const_double_vector  x, complex_vector y) const noexcept = 0;
        /**
         * @brief Discrete forward transformation (without normalization)
         */
        virtual void fwd0(const_double_vector x, complex_vector y) const noexcept = 0;
        /**
         * @brief Discrete forward transformation (unitary transformation, with normalization by sqrt(1/N))
         */
        virtual void fwdu(const_double_vector x, complex_vector y) const noexcept = 0;
        /**
         * @brief Discrete forward transformation (with normalization by 1/N)
         */
        virtual void fwdn(const_double_vector x, complex_vector y) const noexcept = 0;
        /**
         * @brief Discrete inverse transformation (without normalization)
         */
        virtual void inv(complex_vector  x, double_vector y) const noexcept = 0;
        /**
         * @brief Discrete inverse transformation (without normalization)
         */
        virtual void inv0(complex_vector x, double_vector y) const noexcept = 0;
        /**
         * @brief Discrete inverse transformation (unitary transformation)
         */
        virtual void invu(complex_vector x, double_vector y) const noexcept = 0;
        /**
         * @brief Discrete inverse transformation (with normalization by 1/N)
         */
        virtual void invn(complex_vector x, double_vector y) const noexcept = 0;
    };

    /**
     * @brief Real In-Place Discrete Cosine Transform (DCT-II)
     */
    class RealDCT
    {
    public:
        virtual ~RealDCT() = default;

        /**
         * @brief Setup the sequence length of the algorithm (up to 2^30)
         */
        virtual void setup(int n) = 0;

        /**
         * @brief Discrete forward transformation (with normalization by 1/N)
         */
        virtual void fwd(double_vector  x) const noexcept = 0;
        /**
         * @brief Discrete forward transformation (without normalization)
         */
        virtual void fwd0(double_vector x) const noexcept = 0;
        /**
         * @brief Discrete forward transformation (with normalization by 1/N)
         */
        virtual void fwdn(double_vector x) const noexcept = 0;
        /**
         * @brief Discrete inverse transformation (without normalization)
         */
        virtual void inv(double_vector  x) const noexcept = 0;
        /**
         * @brief Discrete inverse transformation (without normalization)
         */
        virtual void inv0(double_vector x) const noexcept = 0;
        /**
         * @brief Discrete inverse transformation (with normalization by 1/N)
         */
        virtual void invn(double_vector x) const noexcept = 0;
    };

    ComplexFFTPtr createComplexFFT(int n);

    RealFFTPtr createRealFFT(int n);

    RealDCTPtr createDCT(int n);

    ComplexFFTPtr createBluesteinFFT(int n);

    /******************************************************************************
    *  Build Info
    ******************************************************************************/
    bool builtWithSSE();
    bool builtWithSSE_MATH();
    bool builtWithSSE2();
    bool builtWithSSE2_MATH();
    bool builtWithSSE3();
    bool builtWithSSSE3();
    bool builtWithSSE4_1();
    bool builtWithSSE4_2();
    bool builtWithAVX();
    bool builtWithAVX2();
    bool builtWithAVX512BW();
    bool builtWithAVX512CD();
    bool builtWithAVX512DQ();
    bool builtWithAVX512F();
    bool builtWithAVX512VL();

    bool builtWithOpenMP();
    int getOpenMPMaxThreads();
}

#endif // otfft_h
