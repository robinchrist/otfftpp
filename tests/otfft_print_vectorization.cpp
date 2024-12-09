#include "otfftpp/otfft.h"
#include <iostream>

int main() {
    if(OTFFT::builtWithOpenMP()) {
        std::cout << "  otfftpp: _OPENMP " << OTFFT::builtWithSSE() << "\n";
        std::cout << "  otfftpp: omp_get_max_threads " << OTFFT::getOpenMPMaxThreads() << "\n";
    }

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