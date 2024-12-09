#include "otfftpp/otfft.h"
#include <memory>
#include <iostream>

int main() {
    auto fft = OTFFT::createBluesteinFFT(13);

    #ifdef NDEBUG
    std::cout << "otfftpp: Release!\n";
    #else
    std::cout << "otfftpp: Debug!\n";
    #endif

    // ARCHITECTURES
    #ifdef _M_X64
    std::cout << "  otfftpp: _M_X64 defined\n";
    #endif

    #ifdef _M_IX86
    std::cout << "  otfftpp: _M_IX86 defined\n";
    #endif

    #ifdef _M_ARM64
    std::cout << "  otfftpp: _M_ARM64 defined\n";
    #endif

    #if __i386__
    std::cout << "  otfftpp: __i386__ defined\n";
    #endif

    #if __x86_64__
    std::cout << "  otfftpp: __x86_64__ defined\n";
    #endif

    #if __aarch64__
    std::cout << "  otfftpp: __aarch64__ defined\n";
    #endif

    // Libstdc++
    #if defined _GLIBCXX_USE_CXX11_ABI
    std::cout << "  otfftpp: _GLIBCXX_USE_CXX11_ABI "<< _GLIBCXX_USE_CXX11_ABI << "\n";
    #endif

    // MSVC runtime
    #if defined(_DEBUG)
        #if defined(_MT) && defined(_DLL)
        std::cout << "  otfftpp: MSVC runtime: MultiThreadedDebugDLL\n";
        #elif defined(_MT)
        std::cout << "  otfftpp: MSVC runtime: MultiThreadedDebug\n";
        #endif
    #else
        #if defined(_MT) && defined(_DLL)
        std::cout << "  otfftpp: MSVC runtime: MultiThreadedDLL\n";
        #elif defined(_MT)
        std::cout << "  otfftpp: MSVC runtime: MultiThreaded\n";
        #endif
    #endif

    // COMPILER VERSIONS
    #if _MSC_VER
    std::cout << "  otfftpp: _MSC_VER " << _MSC_VER<< "\n";
    #endif

    #if _MSVC_LANG
    std::cout << "  otfftpp: _MSVC_LANG " << _MSVC_LANG<< "\n";
    #endif

    #if __cplusplus
    std::cout << "  otfftpp: __cplusplus " << __cplusplus<< "\n";
    #endif

    #if __INTEL_COMPILER
    std::cout << "  otfftpp: __INTEL_COMPILER " << __INTEL_COMPILER<< "\n";
    #endif

    #if __GNUC__
    std::cout << "  otfftpp: __GNUC__ " << __GNUC__<< "\n";
    #endif

    #if __GNUC_MINOR__
    std::cout << "  otfftpp: __GNUC_MINOR__ " << __GNUC_MINOR__<< "\n";
    #endif

    #if __clang_major__
    std::cout << "  otfftpp: __clang_major__ " << __clang_major__<< "\n";
    #endif

    #if __clang_minor__
    std::cout << "  otfftpp: __clang_minor__ " << __clang_minor__<< "\n";
    #endif

    #if __apple_build_version__
    std::cout << "  otfftpp: __apple_build_version__ " << __apple_build_version__<< "\n";
    #endif

    // SUBSYSTEMS

    #if __MSYS__
    std::cout << "  otfftpp: __MSYS__ " << __MSYS__<< "\n";
    #endif

    #if __MINGW32__
    std::cout << "  otfftpp: __MINGW32__ " << __MINGW32__<< "\n";
    #endif

    #if __MINGW64__
    std::cout << "  otfftpp: __MINGW64__ " << __MINGW64__<< "\n";
    #endif

    #if __CYGWIN__
    std::cout << "  otfftpp: __CYGWIN__ " << __CYGWIN__<< "\n";
    #endif

    if(OTFFT::builtWithOpenMP()) {
        std::cout << "  otfftpp: _OPENMP " << OTFFT::builtWithSSE() << "\n";
        std::cout << "  otfftpp: omp_get_max_threads " << OTFFT::getOpenMPMaxThreads() << "\n";
    }

    // INSTRUCTION SETS

    if(OTFFT::builtWithSSE()){
        std::cout << "  otfftpp: __SSE__ " << OTFFT::builtWithSSE() << "\n";
    }

    if(OTFFT::builtWithSSE_MATH()){
        std::cout << "  otfftpp: __SSE_MATH__ " << OTFFT::builtWithSSE_MATH() << "\n";
    }

    if(OTFFT::builtWithSSE2()){
        std::cout << "  otfftpp: __SSE2__ " << OTFFT::builtWithSSE2() << "\n";
    }

    if(OTFFT::builtWithSSE2_MATH()){
        std::cout << "  otfftpp: __SSE2_MATH__ " << OTFFT::builtWithSSE2_MATH() << "\n";
    }

    if(OTFFT::builtWithSSE3()){
        std::cout << "  otfftpp: __SSE3__ " << OTFFT::builtWithSSE3() << "\n";
    }

    if(OTFFT::builtWithSSSE3()){
        std::cout << "  otfftpp: __SSSE3__ " << OTFFT::builtWithSSSE3() << "\n";
    }

    if(OTFFT::builtWithSSE4_1()){
        std::cout << "  otfftpp: __SSE4_1__ " << OTFFT::builtWithSSE4_1() << "\n";
    }

    if(OTFFT::builtWithSSE4_2()){
        std::cout << "  otfftpp: __SSE4_2__ " << OTFFT::builtWithSSE4_2() << "\n";
    }

    if(OTFFT::builtWithAVX()){
        std::cout << "  otfftpp: __AVX__ " << OTFFT::builtWithAVX() << "\n";
    }

    if(OTFFT::builtWithAVX2()){
        std::cout << "  otfftpp: __AVX2__ " << OTFFT::builtWithAVX2() << "\n";
    }

    if(OTFFT::builtWithAVX512BW()){
        std::cout << "  otfftpp: __AVX512BW__ " << OTFFT::builtWithAVX512BW() << "\n";
    }

    if(OTFFT::builtWithAVX512CD()){
        std::cout << "  otfftpp: __AVX512CD__ " << OTFFT::builtWithAVX512CD() << "\n";
    }

    if(OTFFT::builtWithAVX512DQ()){
        std::cout << "  otfftpp: __AVX512DQ__ " << OTFFT::builtWithAVX512DQ() << "\n";
    }

    if(OTFFT::builtWithAVX512F()){
        std::cout << "  otfftpp: __AVX512F__ " << OTFFT::builtWithAVX512F() << "\n";
    }

    if(OTFFT::builtWithAVX512VL()){
        std::cout << "  otfftpp: __AVX512VL__ " << OTFFT::builtWithAVX512VL() << "\n";
}

    
}
