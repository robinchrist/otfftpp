/******************************************************************************
*  OTFFT Implementation Version 11.4xv
*
*  Copyright (c) 2015 OK Ojisan(Takuya OKAHISA)
*  Released under the MIT license
*  http://opensource.org/licenses/mit-license.php
******************************************************************************/

#include "otfftpp/otfft.h"
#include "otfftpp/detail/otfft_misc.h"
#include "otfftpp/detail/otfft_platform.h"

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
#else
#include <unistd.h>
#include <sys/types.h>
#include <sys/param.h>
#include <sys/time.h>
#if defined(PLAT_LINUX)
#elif defined(PLAT_MACOS)
#include <sys/sysctl.h>
#endif
#endif

#if _OPENMP
#include <omp.h>
#endif


#include "otfftpp/detail/otfft_if.h"
#include "otfftpp/detail/otfft_base.h"


namespace OTFFT
{
    ComplexFFTPtr createComplexFFT(int n)
    {
        return ComplexFFTPtr(new FFT(n));
    }

    RealFFTPtr createRealFFT(int n)
    {
        return RealFFTPtr(new RFFT(n));
    }

    RealDCTPtr createDCT(int n)
    {
        return RealDCTPtr(new DCT(n));
    }

    ComplexFFTPtr createBluesteinFFT(int n)
    {
        return ComplexFFTPtr(new Bluestein(n));
    }

    bool builtWithSSE() {
        #if __SSE__
        return true;
        #endif
        return false;
    }
    bool builtWithSSE_MATH() {
        #if __SSE_MATH__
        return true;
        #endif
        return false;
    }
    bool builtWithSSE2() {
        #if __SSE2__
        return true;
        #endif
        return false;
    }
    bool builtWithSSE2_MATH() {
        #if __SSE2_MATH__
        return true;
        #endif
        return false;
    }
    bool builtWithSSE3() {
        #if __SSE3__
        return true;
        #endif
        return false;
    }
    bool builtWithSSSE3() {
        #if __SSSE3__
        return true;
        #endif
        return false;
    }
    bool builtWithSSE4_1() {
        #if __SSE4_1__
        return true;
        #endif
        return false;
    }
    bool builtWithSSE4_2() {
        #if __SSE4_2__
        return true;
        #endif
        return false;
    }
    bool builtWithAVX() {
        #if __AVX__
        return true;
        #endif
        return false;
    }
    bool builtWithAVX2() {
        #if __AVX2__
        return true;
        #endif
        return false;
    }
    bool builtWithAVX512BW() {
        #if __AVX512BW__
        return true;
        #endif
        return false;
    }
    bool builtWithAVX512CD() {
        #if __AVX512CD__
        return true;
        #endif
        return false;
    }
    bool builtWithAVX512DQ() {
        #if __AVX512DQ__
        return true;
        #endif
        return false;
    }
    bool builtWithAVX512F() {
        #if __AVX512F__
        return true;
        #endif
        return false;
    }
    bool builtWithAVX512VL() {
        #if __AVX512VL__
        return true;
        #endif
        return false;
    }

    bool builtWithOpenMP() {
        #if _OPENMP
        return true;
        #endif
        return false;
    }

    int getOpenMPMaxThreads() {
        #if _OPENMP
        return omp_get_max_threads();
        #endif
        return -1;
    }
}
