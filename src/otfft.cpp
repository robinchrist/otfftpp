/******************************************************************************
*  OTFFT Implementation Version 11.4xv
*
*  Copyright (c) 2015 OK Ojisan(Takuya OKAHISA)
*  Released under the MIT license
*  http://opensource.org/licenses/mit-license.php
******************************************************************************/

#include "otfft_config.h"
#include "otfft.h"
#include "otfft_misc.h"
#include "otfft_platform.h"

#include <thread>
#include <array>
#include <stdexcept>
#include <cstdint>
#include <cassert>

#include <stdlib.h>
#include <string.h>

#if defined(PLAT_MINGW)
#include <cpuid.h>
#include <windows.h>
#elif defined(PLAT_WINDOWS)
#include <windows.h>
#include <limits.h>
#include <intrin.h>
#include <immintrin.h>
#else
#include <unistd.h>
#include <sys/types.h>
#include <sys/param.h>
#include <sys/time.h>
#if defined(PLAT_LINUX)
#include <x86intrin.h>
#elif defined(PLAT_MACOS)
#include <sys/sysctl.h>
#endif
#endif


#ifdef OTFFT_WITH_AVX2
#include "otfft_avx2.h"
#endif

namespace OTFFT
{
    void unique_ptr_delete::operator()(ComplexFFT *raw_pointer)
    {
        return OTFFT_AVX2::unique_ptr_deleter(raw_pointer);
    }

    void unique_ptr_delete::operator()(RealFFT *raw_pointer)
    {
        return OTFFT_AVX2::unique_ptr_deleter(raw_pointer);
    }

    void unique_ptr_delete::operator()(RealDCT *raw_pointer)
    {
        return OTFFT_AVX2::unique_ptr_deleter(raw_pointer);
    }

    namespace Factory
    {
        ComplexFFTPtr createComplexFFT(int n, OptimizationType t)
        {
            return ComplexFFTPtr(OTFFT_AVX2::Factory::createComplexFFT(n));
        }

        RealFFTPtr createRealFFT(int n, OptimizationType t)
        {
            return RealFFTPtr(OTFFT_AVX2::Factory::createRealFFT(n));
        }

        RealDCTPtr createDCT(int n, OptimizationType t)
        {
            return RealDCTPtr(OTFFT_AVX2::Factory::createDCT(n));
        }

        ComplexFFTPtr createBluesteinFFT(int n, OptimizationType t)
        {
            return ComplexFFTPtr(OTFFT_AVX2::Factory::createBluesteinFFT(n));
        }
    }
}
