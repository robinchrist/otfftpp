// Copyright (C) DEWETRON GmbH 2017

#pragma once

#include <memory>

namespace OTFFT
{
    struct complex_t;

    class ComplexFFT;
    class RealFFT;
    class RealDCT;

    using ComplexFFTPtr = std::unique_ptr<ComplexFFT>;
    using RealFFTPtr = std::unique_ptr<RealFFT>;
    using RealDCTPtr = std::unique_ptr<RealDCT>;
}
