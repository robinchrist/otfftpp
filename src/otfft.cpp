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


#include "otfft_if.h"
#include "otfft_base.h"


namespace OTFFT
{
    namespace Factory
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
    }
}
