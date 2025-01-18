/******************************************************************************
*  FFT Miscellaneous Routines Version 11.4xv
*
*  Copyright (c) 2019 OK Ojisan(Takuya OKAHISA)
*  Released under the MIT license
*  http://opensource.org/licenses/mit-license.php
******************************************************************************/
// Copyright (c) 2017 to the present, DEWETRON GmbH

#ifndef otfft_misc_h
#define otfft_misc_h

//=============================================================================
// Customization Options
//=============================================================================

#define USE_INTRINSIC 1
//#define DO_SINGLE_THREAD 1
//#define USE_UNALIGNED_MEMORY 1

//=============================================================================

#include "otfftpp/otfft_complex.h"
#include <complex>
#include <cmath>
#include <new>

#include "simde/x86/sse2.h"
#include "simde/x86/sse3.h"
#include "simde/x86/avx.h"
#include "simde/x86/avx512.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

#ifndef M_SQRT2
#define M_SQRT2 1.41421356237309504876378807303183294
#endif

#ifndef M_SQRT1_2
#define M_SQRT1_2 0.707106781186547524400844362104849039
#endif

namespace OTFFT_MISC {
    using namespace OTFFT;

constexpr double H1X =  0.923879532511286762010323247995557949;
constexpr double H1Y = -0.382683432365089757574419179753100195;

enum scaling_mode { scale_1 = 0, scale_unitary, scale_length };

static inline complex_t v8x(const complex_t& z) noexcept force_inline;
static inline complex_t v8x(const complex_t& z) noexcept
{
    return complex_t(M_SQRT1_2*(z.Re-z.Im), M_SQRT1_2*(z.Re+z.Im));
}
static inline complex_t w8x(const complex_t& z) noexcept force_inline;
static inline complex_t w8x(const complex_t& z) noexcept
{
    return complex_t(M_SQRT1_2*(z.Re+z.Im), M_SQRT1_2*(z.Im-z.Re));
}

} // namespace OTFFT_MISC


//=============================================================================
// constexpr sqrt
//=============================================================================

namespace OTFFT_MISC {

constexpr double sqrt_aux(double a, double x, double y)
{
    return x == y ? x : sqrt_aux(a, (x + a/x)/2, x);
}

constexpr double mysqrt(double x) { return sqrt_aux(x, x/2, x); }

constexpr int mylog2(int N)
{
    return N <= 1 ? 0 : 1 + mylog2(N/2);
}

template <int N, int s>
static complex_t modq(const_complex_vector W, const int p) noexcept
{
    constexpr int Nq = N/4;
    constexpr int log_Nq = mylog2(Nq);
    const int sp = s*p;
    const int q = sp >> log_Nq;
    const int r = sp & (Nq-1);
    const complex_t z = W[r];
    switch (q & 3) {
        case 0: return z;
        case 1: return mjx(z);
        case 2: return neg(z);
        case 3: return jx(z);
    }
    return complex_t();
}

} // namespace OTFFT_MISC

//=============================================================================
// FFT Weight Initialize Routine
//=============================================================================

namespace OTFFT_MISC {

#ifdef DO_SINGLE_THREAD
constexpr int OMP_THRESHOLD_W = 1<<30;
#else
//constexpr int OMP_THRESHOLD_W = 1<<16;
constexpr int OMP_THRESHOLD_W = 1<<15;
#endif

    static void init_W(const int N, complex_vector W) noexcept
    {
        const double theta0 = 2*M_PI/N;
        const int Nh = N/2;
        const int Nq = N/4;
        const int Ne = N/8;
        const int Nd = N - Nq;
        if (N < 1) {}
        else if (N < 2) { W[0] = W[1] = 1; }
        else if (N < 4) { W[0] = W[2] = 1; W[1] = -1; }
        else if (N < 8) {
            W[0] = complex_t( 1,  0);
            W[1] = complex_t( 0, -1);
            W[2] = complex_t(-1,  0);
            W[3] = complex_t( 0,  1);
            W[4] = complex_t( 1,  0);
        }
        else if (N < OMP_THRESHOLD_W) for (int p = 0; p <= Ne; p++) {
            const double theta = p * theta0;
            const double c =  cos(theta);
            const double s = -sin(theta);
            W[p]    = complex_t( c,  s);
            W[Nq-p] = complex_t(-s, -c);
            W[Nq+p] = complex_t( s, -c);
            W[Nh-p] = complex_t(-c,  s);
            W[Nh+p] = complex_t(-c, -s);
            W[Nd-p] = complex_t( s,  c);
            W[Nd+p] = complex_t(-s,  c);
            W[N-p]  = complex_t( c, -s);
        }
        else {
        #pragma omp parallel for schedule(static)
        for (int p = 0; p <= Ne; p++) {
            const double theta = p * theta0;
            const double c =  cos(theta);
            const double s = -sin(theta);
            W[p]    = complex_t( c,  s);
            W[Nq-p] = complex_t(-s, -c);
            W[Nq+p] = complex_t( s, -c);
            W[Nh-p] = complex_t(-c,  s);
            W[Nh+p] = complex_t(-c, -s);
            W[Nd-p] = complex_t( s,  c);
            W[Nd+p] = complex_t(-s,  c);
            W[N-p]  = complex_t( c, -s);
        }
        }
    }

    static void init_Wr1(const int r, const int N, complex_vector W) noexcept
    {
       if (N < r) return;
       const int Nr = N/r;
       const double theta = -2*M_PI/N;
       if (N < OMP_THRESHOLD_W) {
           for (int p = 0; p < Nr; p++) {
               W[p] = expj(theta * p);
           }
       }
       else {
#pragma omp parallel for schedule(static)
           for (int p = 0; p < Nr; p++) {
               W[p] = expj(theta * p);
           }
       }
    }

    static void init_Wt(const int r, const int N, complex_vector W) noexcept
    {
       if (N < r) return;
       const int Nr = N/r;
       const double theta = -2*M_PI/N;
       if (N < OMP_THRESHOLD_W) {
           for (int p = 0; p < Nr; p++) {
               for (int k = 1; k < r; k++) {
                   W[p + (k-1)*Nr] = W[N + r*p + k] = expj(theta * k*p);
               }
           }
       }
       else {
#pragma omp parallel for schedule(static)
           for (int p = 0; p < Nr; p++) {
               for (int k = 1; k < r; k++) {
                   W[p + (k-1)*Nr] = W[N + r*p + k] = expj(theta * k*p);
               }
           }
       }
    }

    template <int r, int N, int k>
    const complex_t* twid(const_complex_vector W, int p)
    {
        constexpr int Nr = N/r;
        constexpr int d = (k-1)*Nr;
        return W + p + d;
    }
    
    template <int r, int N, int k>
    const complex_t* twidT(const_complex_vector W, int p)
    {
        constexpr int d = N + k;
        return W + r*p + d;
    }

    static void speedup_magic(const int N = 1 << 18) noexcept
    {
        const double theta0 = 2*M_PI/N;
        volatile double sum = 0;
        for (int p = 0; p < N; p++) {
            sum += cos(p * theta0);
        }
    }

} // namespace OTFFT_MISC

namespace OTFFT_MISC {

    using xmm = simde__m128d;

    static inline xmm cmplx(const double& x, const double& y) noexcept force_inline;
    static inline xmm cmplx(const double& x, const double& y) noexcept
    {
        return simde_mm_setr_pd(x, y);
    }

    static inline xmm getpz(const complex_t& z) noexcept force_inline;
    static inline xmm getpz(const complex_t& z) noexcept
    {
#ifdef USE_UNALIGNED_MEMORY
        return simde_mm_loadu_pd(&z.Re);
#else
        return simde_mm_load_pd(&z.Re);
#endif
    }
    static inline xmm getpz(const_double_vector x) noexcept force_inline;
    static inline xmm getpz(const_double_vector x) noexcept
    {
#ifdef USE_UNALIGNED_MEMORY
        return simde_mm_loadu_pd(x);
#else
        return simde_mm_load_pd(x);
#endif
    }

    static inline void setpz(complex_t& z, const xmm x) noexcept force_inline3;
    static inline void setpz(complex_t& z, const xmm x) noexcept
    {
#ifdef USE_UNALIGNED_MEMORY
        simde_mm_storeu_pd(&z.Re, x);
#else
        simde_mm_store_pd(&z.Re, x);
#endif
    }
    static inline void setpz(double_vector x, const xmm z) noexcept force_inline3;
    static inline void setpz(double_vector x, const xmm z) noexcept
    {
#ifdef USE_UNALIGNED_MEMORY
        simde_mm_storeu_pd(x, z);
#else
        simde_mm_store_pd(x, z);
#endif
    }

    static inline void swappz(complex_t& x, complex_t& y) noexcept force_inline3;
    static inline void swappz(complex_t& x, complex_t& y) noexcept
    {
        const xmm z = getpz(x); setpz(x, getpz(y)); setpz(y, z);
    }

    static inline xmm cnjpz(const xmm xy) noexcept force_inline;
    static inline xmm cnjpz(const xmm xy) noexcept
    {
        constexpr xmm zm = { 0.0, -0.0 };
        return simde_mm_xor_pd(zm, xy);
    }
    static inline xmm jxpz(const xmm xy) noexcept force_inline;
    static inline xmm jxpz(const xmm xy) noexcept
    {
        const xmm xmy = cnjpz(xy);
        return simde_mm_shuffle_pd(xmy, xmy, 1);
    }
    static inline xmm negpz(const xmm xy) noexcept force_inline;
    static inline xmm negpz(const xmm xy) noexcept
    {
        constexpr xmm mm = { -0.0, -0.0 };
        return simde_mm_xor_pd(mm, xy);
    }
    static inline xmm mjxpz(const xmm xy) noexcept force_inline;
    static inline xmm mjxpz(const xmm xy) noexcept
    {
        const xmm yx = simde_mm_shuffle_pd(xy, xy, 1);
        return cnjpz(yx);
    }

    static inline xmm addpz(const xmm a, const xmm b) noexcept force_inline;
    static inline xmm addpz(const xmm a, const xmm b) noexcept
    {
        return simde_mm_add_pd(a, b);
    }
    static inline xmm subpz(const xmm a, const xmm b) noexcept force_inline;
    static inline xmm subpz(const xmm a, const xmm b) noexcept
    {
        return simde_mm_sub_pd(a, b);
    }
    static inline xmm mulpd(const xmm a, const xmm b) noexcept force_inline;
    static inline xmm mulpd(const xmm a, const xmm b) noexcept
    {
        return simde_mm_mul_pd(a, b);
    }
    static inline xmm divpd(const xmm a, const xmm b) noexcept force_inline;
    static inline xmm divpd(const xmm a, const xmm b) noexcept
    {
        return simde_mm_div_pd(a, b);
    }

    template <int N, int mode> static inline xmm scalepz(const xmm z) noexcept force_inline;
    template <int N, int mode> static inline xmm scalepz(const xmm z) noexcept
    {
        constexpr double scale =
            mode == scale_1       ? 1.0           :
            mode == scale_unitary ? 1.0/mysqrt(N) :
            mode == scale_length  ? 1.0/N         : 0.0;
        constexpr xmm sv = { scale, scale };
        return mode == scale_1 ? z : mulpd(sv, z);
    }

template <int N, int s>
static xmm modqpz(const_complex_vector W, const int p) noexcept
{
    constexpr int Nq = N/4;
    constexpr int log_Nq = mylog2(Nq);
    const int sp = s*p;
    const int q = sp >> log_Nq;
    const int r = sp & (Nq-1);
    const xmm x = getpz(W[r]);
    switch (q & 3) {
        case 0: return x;
        case 1: return mjxpz(x);
        case 2: return negpz(x);
        case 3: return jxpz(x);
    }
    return xmm();
}

} // namespace OTFFT_MISC

//-----------------------------------------------------------------------------
// SSE3
//-----------------------------------------------------------------------------

namespace OTFFT_MISC {

    static inline xmm haddpz(const xmm ab, const xmm xy) noexcept force_inline;
    static inline xmm haddpz(const xmm ab, const xmm xy) noexcept
    {
        return simde_mm_hadd_pd(ab, xy); // (a + b, x + y)
    }

    static inline xmm mulpz(const xmm ab, const xmm xy) noexcept force_inline;
    static inline xmm mulpz(const xmm ab, const xmm xy) noexcept
    {
        const xmm aa = simde_mm_unpacklo_pd(ab, ab);
        const xmm bb = simde_mm_unpackhi_pd(ab, ab);
        const xmm yx = simde_mm_shuffle_pd(xy, xy, 1);
        return simde_mm_fmaddsub_pd(aa, xy, simde_mm_mul_pd(bb, yx));
    }

    static inline xmm divpz(const xmm ab, const xmm xy) noexcept force_inline;
    static inline xmm divpz(const xmm ab, const xmm xy) noexcept
    {
        const xmm x2y2 = simde_mm_mul_pd(xy, xy);
        const xmm r2r2 = simde_mm_hadd_pd(x2y2, x2y2);
        return simde_mm_div_pd(mulpz(ab, cnjpz(xy)), r2r2);
    }

    static inline xmm v8xpz(const xmm xy) noexcept force_inline;
    static inline xmm v8xpz(const xmm xy) noexcept
    {
        constexpr xmm rr = { M_SQRT1_2, M_SQRT1_2 };
        const xmm yx = simde_mm_shuffle_pd(xy, xy, 1);
        return simde_mm_mul_pd(rr, simde_mm_addsub_pd(xy, yx));
    }

} // namespace OTFFT_MISC

namespace OTFFT_MISC {

    static inline xmm w8xpz(const xmm xy) noexcept force_inline;
    static inline xmm w8xpz(const xmm xy) noexcept
    {
        constexpr xmm rr = { M_SQRT1_2, M_SQRT1_2 };
        const xmm ymx = cnjpz(simde_mm_shuffle_pd(xy, xy, 1));
        return simde_mm_mul_pd(rr, simde_mm_add_pd(xy, ymx));
    }

    static inline xmm h1xpz(const xmm xy) noexcept force_inline;
    static inline xmm h1xpz(const xmm xy) noexcept
    {
        constexpr xmm h1 = { H1X, H1Y };
        return mulpz(h1, xy);
    }

    static inline xmm h3xpz(const xmm xy) noexcept force_inline;
    static inline xmm h3xpz(const xmm xy) noexcept
    {
        constexpr xmm h3 = { -H1Y, -H1X };
        return mulpz(h3, xy);
    }

    static inline xmm hfxpz(const xmm xy) noexcept force_inline;
    static inline xmm hfxpz(const xmm xy) noexcept
    {
        constexpr xmm hf = { H1X, -H1Y };
        return mulpz(hf, xy);
    }

    static inline xmm hdxpz(const xmm xy) noexcept force_inline;
    static inline xmm hdxpz(const xmm xy) noexcept
    {
        constexpr xmm hd = { -H1Y, H1X };
        return mulpz(hd, xy);
    }
} // namespace OTFFT_MISC



namespace OTFFT_MISC {

    using ymm = simde__m256d;

    static inline ymm cmplx2(const double a, const double b, const double c, const double d) noexcept force_inline;
    static inline ymm cmplx2(const double a, const double b, const double c, const double d) noexcept
    {
        return simde_mm256_setr_pd(a, b, c, d);
    }

    static inline ymm cmplx2(const complex_t& x, const complex_t& y) noexcept force_inline;
    static inline ymm cmplx2(const complex_t& x, const complex_t& y) noexcept
    {
        return simde_mm256_setr_pd(x.Re, x.Im, y.Re, y.Im);
    }

    static inline ymm getpz2(const_complex_vector z) noexcept force_inline;
    static inline ymm getpz2(const_complex_vector z) noexcept
    {
#ifdef USE_UNALIGNED_MEMORY
        return simde_mm256_loadu_pd(&z->Re);
#else
        return simde_mm256_load_pd(&z->Re);
#endif
    }

    static inline void setpz2(complex_vector z, const ymm x) noexcept force_inline3;
    static inline void setpz2(complex_vector z, const ymm x) noexcept
    {
#ifdef USE_UNALIGNED_MEMORY
        simde_mm256_storeu_pd(&z->Re, x);
#else
        simde_mm256_store_pd(&z->Re, x);
#endif
    }

    static inline ymm cnjpz2(const ymm xy) noexcept force_inline;
    static inline ymm cnjpz2(const ymm xy) noexcept
    {
        constexpr ymm zm = { 0.0, -0.0, 0.0, -0.0 };
        return simde_mm256_xor_pd(zm, xy);
    }
    static inline ymm jxpz2(const ymm xy) noexcept force_inline;
    static inline ymm jxpz2(const ymm xy) noexcept
    {
        const ymm xmy = cnjpz2(xy);
        return simde_mm256_shuffle_pd(xmy, xmy, 5);
    }
    static inline ymm negpz2(const ymm xy) noexcept force_inline;
    static inline ymm negpz2(const ymm xy) noexcept
    {
        constexpr ymm mm = { -0.0, -0.0, -0.0, -0.0 };
        return simde_mm256_xor_pd(mm, xy);
    }

    static inline ymm addpz2(const ymm a, const ymm b) noexcept force_inline;
    static inline ymm addpz2(const ymm a, const ymm b) noexcept
    {
        return simde_mm256_add_pd(a, b);
    }
    static inline ymm subpz2(const ymm a, const ymm b) noexcept force_inline;
    static inline ymm subpz2(const ymm a, const ymm b) noexcept
    {
        return simde_mm256_sub_pd(a, b);
    }
    static inline ymm mulpd2(const ymm a, const ymm b) noexcept force_inline;
    static inline ymm mulpd2(const ymm a, const ymm b) noexcept
    {
        return simde_mm256_mul_pd(a, b);
    }
    static inline ymm divpd2(const ymm a, const ymm b) noexcept force_inline;
    static inline ymm divpd2(const ymm a, const ymm b) noexcept
    {
        return simde_mm256_div_pd(a, b);
    }

    static inline ymm mulpz2(const ymm ab, const ymm xy) noexcept force_inline;
    static inline ymm mulpz2(const ymm ab, const ymm xy) noexcept
    {
        const ymm aa = simde_mm256_unpacklo_pd(ab, ab);
        const ymm bb = simde_mm256_unpackhi_pd(ab, ab);
        const ymm yx = simde_mm256_shuffle_pd(xy, xy, 5);

        return simde_mm256_fmaddsub_pd(aa, xy, simde_mm256_mul_pd(bb, yx));
    }

    static inline ymm divpz2(const ymm ab, const ymm xy) noexcept force_inline;
    static inline ymm divpz2(const ymm ab, const ymm xy) noexcept
    {
        const ymm x2y2 = simde_mm256_mul_pd(xy, xy);
        const ymm r2r2 = simde_mm256_hadd_pd(x2y2, x2y2);
        return simde_mm256_div_pd(mulpz2(ab, cnjpz2(xy)), r2r2);
    }

    template <int N, int mode> static inline ymm scalepz2(const ymm z) noexcept force_inline;
    template <int N, int mode> static inline ymm scalepz2(const ymm z) noexcept
    {
        constexpr double scale =
            mode == scale_1       ? 1.0           :
            mode == scale_unitary ? 1.0/mysqrt(N) :
            mode == scale_length  ? 1.0/N         : 0.0;
        constexpr ymm sv = { scale, scale, scale, scale };
        return mode == scale_1 ? z : mulpd2(sv, z);
    }

    static inline ymm v8xpz2(const ymm xy) noexcept force_inline;
    static inline ymm v8xpz2(const ymm xy) noexcept
    {
        constexpr ymm rr = { M_SQRT1_2, M_SQRT1_2, M_SQRT1_2, M_SQRT1_2 };
        const ymm yx = simde_mm256_shuffle_pd(xy, xy, 5);
        return simde_mm256_mul_pd(rr, simde_mm256_addsub_pd(xy, yx));
    }

    static inline ymm w8xpz2(const ymm xy) noexcept force_inline;
    static inline ymm w8xpz2(const ymm xy) noexcept
    {
        constexpr ymm rr = { M_SQRT1_2, M_SQRT1_2, M_SQRT1_2, M_SQRT1_2 };
        const ymm ymx = cnjpz2(simde_mm256_shuffle_pd(xy, xy, 5));
        return simde_mm256_mul_pd(rr, simde_mm256_add_pd(xy, ymx));
    }

    static inline ymm h1xpz2(const ymm xy) noexcept force_inline;
    static inline ymm h1xpz2(const ymm xy) noexcept
    {
        constexpr ymm h1 = { H1X, H1Y, H1X, H1Y };
        return mulpz2(h1, xy);
    }

    static inline ymm h3xpz2(const ymm xy) noexcept force_inline;
    static inline ymm h3xpz2(const ymm xy) noexcept
    {
        constexpr ymm h3 = { -H1Y, -H1X, -H1Y, -H1X };
        return mulpz2(h3, xy);
    }

    static inline ymm hfxpz2(const ymm xy) noexcept force_inline;
    static inline ymm hfxpz2(const ymm xy) noexcept
    {
        constexpr ymm hf = { H1X, -H1Y, H1X, -H1Y };
        return mulpz2(hf, xy);
    }

    static inline ymm hdxpz2(const ymm xy) noexcept force_inline;
    static inline ymm hdxpz2(const ymm xy) noexcept
    {
        constexpr ymm hd = { -H1Y, H1X, -H1Y, H1X };
        return mulpz2(hd, xy);
    }

    static inline ymm duppz2(const xmm x) noexcept force_inline;
    static inline ymm duppz2(const xmm x) noexcept
    {
        return simde_mm256_broadcast_pd(&x);
    }

    static inline ymm duppz3(const complex_t& z) noexcept force_inline;
    static inline ymm duppz3(const complex_t& z) noexcept
    {
        return simde_mm256_broadcast_pd(reinterpret_cast<const xmm *>(&z));
    }

    static inline ymm cat(const xmm a, const xmm b) noexcept force_inline;
    static inline ymm cat(const xmm a, const xmm b) noexcept
    {
        const ymm ax = simde_mm256_castpd128_pd256(a);
        return simde_mm256_insertf128_pd(ax, b, 1);
    }

    static inline ymm catlo(const ymm ax, const ymm by) noexcept force_inline;
    static inline ymm catlo(const ymm ax, const ymm by) noexcept
    {
        return simde_mm256_permute2f128_pd(ax, by, 0x20); // == ab
    }

    static inline ymm cathi(const ymm ax, const ymm by) noexcept force_inline;
    static inline ymm cathi(const ymm ax, const ymm by) noexcept
    {
        return simde_mm256_permute2f128_pd(ax, by, 0x31); // == xy
    }

    static inline ymm swaplohi(const ymm ab) noexcept force_inline;
    static inline ymm swaplohi(const ymm ab) noexcept
    {
        return simde_mm256_permute2f128_pd(ab, ab, 0x01); // == ba
    }

    template <int s> static inline ymm getwp2(const_complex_vector W, const int p) noexcept force_inline;
    template <int s> static inline ymm getwp2(const_complex_vector W, const int p) noexcept
    {
        const int sp = s*p;
        return cmplx2(W[sp], W[sp+s]);
    }

    template <int s> static inline ymm cnj_getwp2(const_complex_vector W, const int p) noexcept force_inline;
    template <int s> static inline ymm cnj_getwp2(const_complex_vector W, const int p) noexcept
    {
        const int sp = s*p;
        return cnjpz2(cmplx2(W[sp], W[sp+s]));
    }

    static inline xmm getlo(const ymm a_b) noexcept force_inline;
    static inline xmm getlo(const ymm a_b) noexcept
    {
        return simde_mm256_castpd256_pd128(a_b); // == a
    }
    static inline xmm gethi(const ymm a_b) noexcept force_inline;
    static inline xmm gethi(const ymm a_b) noexcept
    {
        return simde_mm256_extractf128_pd(a_b, 1); // == b
    }

    template <int s> static inline ymm getpz3(const_complex_vector z) noexcept force_inline;
    template <int s> static inline ymm getpz3(const_complex_vector z) noexcept
    {
        return cmplx2(z[0], z[s]);
    }

    template <int s> static inline void setpz3(complex_vector z, const ymm x) noexcept force_inline3;
    template <int s> static inline void setpz3(complex_vector z, const ymm x) noexcept
    {
        setpz(z[0], getlo(x));
        setpz(z[s], gethi(x));
    }

} // namespace OTFFT_MISC

#if defined(__AVX512F__) && defined(__AVX512DQ__) && defined(USE_INTRINSIC)
//=============================================================================
// AVX-512
//=============================================================================

namespace OTFFT_MISC {

using zmm = simde__m512d;

static inline zmm getpz4(const_complex_vector z) noexcept force_inline;
static inline zmm getpz4(const_complex_vector z) noexcept
{
#ifdef USE_UNALIGNED_MEMORY
    return simde_mm512_loadu_pd(&z->Re);
#else
    return simde_mm512_load_pd(&z->Re);
#endif
}

static inline void setpz4(complex_vector a, const zmm z) noexcept force_inline3;
static inline void setpz4(complex_vector a, const zmm z) noexcept
{
#ifdef USE_UNALIGNED_MEMORY
    simde_mm512_storeu_pd(&a->Re, z);
#else
    simde_mm512_store_pd(&a->Re, z);
#endif
}

static inline zmm cnjpz4(const zmm xy) noexcept force_inline;
static inline zmm cnjpz4(const zmm xy) noexcept
{
    constexpr zmm zm = { 0.0, -0.0, 0.0, -0.0, 0.0, -0.0, 0.0, -0.0 };
    return simde_mm512_xor_pd(zm, xy);
}
static inline zmm jxpz4(const zmm xy) noexcept force_inline;
static inline zmm jxpz4(const zmm xy) noexcept
{
    const zmm xmy = cnjpz4(xy);
    return simde_mm512_shuffle_pd(xmy, xmy, 0x55);
}

static inline zmm addpz4(const zmm a, const zmm b) noexcept force_inline;
static inline zmm addpz4(const zmm a, const zmm b) noexcept
{
    return simde_mm512_add_pd(a, b);
}
static inline zmm subpz4(const zmm a, const zmm b) noexcept force_inline;
static inline zmm subpz4(const zmm a, const zmm b) noexcept
{
    return simde_mm512_sub_pd(a, b);
}
static inline zmm mulpd4(const zmm a, const zmm b) noexcept force_inline;
static inline zmm mulpd4(const zmm a, const zmm b) noexcept
{
    return simde_mm512_mul_pd(a, b);
}
static inline zmm divpd4(const zmm a, const zmm b) noexcept force_inline;
static inline zmm divpd4(const zmm a, const zmm b) noexcept
{
    return simde_mm512_div_pd(a, b);
}

static inline zmm mulpz4(const zmm ab, const zmm xy) noexcept force_inline;
static inline zmm mulpz4(const zmm ab, const zmm xy) noexcept
{
    const zmm aa = simde_mm512_unpacklo_pd(ab, ab);
    const zmm bb = simde_mm512_unpackhi_pd(ab, ab);
    const zmm yx = simde_mm512_shuffle_pd(xy, xy, 0x55);
    //TODO: This instruction does not yet exist in SIMDe, an issue has been filed.
    return simde_mm512_fmaddsub_pd(aa, xy, simde_mm512_mul_pd(bb, yx));
}

template <int N, int mode> static inline zmm scalepz4(const zmm z) noexcept force_inline;
template <int N, int mode> static inline zmm scalepz4(const zmm z) noexcept
{
    constexpr double sc =
        mode == scale_1       ? 1.0           :
        mode == scale_unitary ? 1.0/mysqrt(N) :
        mode == scale_length  ? 1.0/N         : 0.0;
    constexpr zmm sv  = { sc, sc, sc, sc, sc, sc, sc, sc };
    return mode == scale_1 ? z : mulpd4(sv, z);
}

static inline zmm v8xpz4(const zmm xy) noexcept force_inline;
static inline zmm v8xpz4(const zmm xy) noexcept
{
    constexpr double r = M_SQRT1_2;
    constexpr zmm rr = { r, r, r, r, r, r, r, r };
    const zmm myx = jxpz4(xy);
    return simde_mm512_mul_pd(rr, simde_mm512_add_pd(xy, myx));
}

static inline zmm w8xpz4(const zmm xy) noexcept force_inline;
static inline zmm w8xpz4(const zmm xy) noexcept
{
    constexpr double r = M_SQRT1_2;
    constexpr zmm rr = { r, r, r, r, r, r, r, r };
    const zmm ymx = cnjpz4(simde_mm512_shuffle_pd(xy, xy, 0x55));
    return simde_mm512_mul_pd(rr, simde_mm512_add_pd(xy, ymx));
}

static inline zmm h1xpz4(const zmm xy) noexcept force_inline;
static inline zmm h1xpz4(const zmm xy) noexcept
{
    constexpr zmm h1 = { H1X, H1Y, H1X, H1Y, H1X, H1Y, H1X, H1Y };
    return mulpz4(h1, xy);
}

static inline zmm h3xpz4(const zmm xy) noexcept force_inline;
static inline zmm h3xpz4(const zmm xy) noexcept
{
    constexpr zmm h3 = { -H1Y, -H1X, -H1Y, -H1X, -H1Y, -H1X, -H1Y, -H1X };
    return mulpz4(h3, xy);
}

static inline zmm hfxpz4(const zmm xy) noexcept force_inline;
static inline zmm hfxpz4(const zmm xy) noexcept
{
    constexpr zmm hf = { H1X, -H1Y, H1X, -H1Y, H1X, -H1Y, H1X, -H1Y };
    return mulpz4(hf, xy);
}

static inline zmm hdxpz4(const zmm xy) noexcept force_inline;
static inline zmm hdxpz4(const zmm xy) noexcept
{
    constexpr zmm hd = { -H1Y, H1X, -H1Y, H1X, -H1Y, H1X, -H1Y, H1X };
    return mulpz4(hd, xy);
}

static inline zmm duppz4(const xmm x) noexcept force_inline;
static inline zmm duppz4(const xmm x) noexcept
{
    return simde_mm512_broadcast_f64x2(x);
}

static inline zmm duppz5(const complex_t& z) noexcept force_inline;
static inline zmm duppz5(const complex_t& z) noexcept
{
    return duppz4(getpz(z));
}

} // namespace OTFFT_MISC

#else
//TODO: Remove as soon as _mm512_fmaddsub_pd is implemented in SIMDe
//=============================================================================
// AVX-512 Emulation
//=============================================================================

namespace OTFFT_MISC {

struct zmm { ymm lo, hi; };

static inline zmm getpz4(const_complex_vector a) noexcept force_inline;
static inline zmm getpz4(const_complex_vector a) noexcept
{
    const zmm z = { getpz2(&a[0]), getpz2(&a[2]) };
    return z;
}

static inline void setpz4(complex_vector a, const zmm& z) noexcept force_inline3;
static inline void setpz4(complex_vector a, const zmm& z) noexcept
{
    setpz2(&a[0], z.lo);
    setpz2(&a[2], z.hi);
}

static inline zmm cnjpz4(const zmm& xy) noexcept force_inline;
static inline zmm cnjpz4(const zmm& xy) noexcept
{
    const zmm z = { cnjpz2(xy.lo), cnjpz2(xy.hi) };
    return z;
}
static inline zmm jxpz4(const zmm& xy) noexcept force_inline;
static inline zmm jxpz4(const zmm& xy) noexcept
{
    const zmm z = { jxpz2(xy.lo), jxpz2(xy.hi) };
    return z;
}

static inline zmm addpz4(const zmm& a, const zmm& b) noexcept force_inline;
static inline zmm addpz4(const zmm& a, const zmm& b) noexcept
{
    const zmm z = { addpz2(a.lo, b.lo), addpz2(a.hi, b.hi) };
    return z;
}
static inline zmm subpz4(const zmm& a, const zmm& b) noexcept force_inline;
static inline zmm subpz4(const zmm& a, const zmm& b) noexcept
{
    const zmm z = { subpz2(a.lo, b.lo), subpz2(a.hi, b.hi) };
    return z;
}
static inline zmm mulpd4(const zmm& a, const zmm& b) noexcept force_inline;
static inline zmm mulpd4(const zmm& a, const zmm& b) noexcept
{
    const zmm z = { mulpd2(a.lo, b.lo), mulpd2(a.hi, b.hi) };
    return z;
}
static inline zmm divpd4(const zmm& a, const zmm& b) noexcept force_inline;
static inline zmm divpd4(const zmm& a, const zmm& b) noexcept
{
    const zmm z = { divpd2(a.lo, b.lo), divpd2(a.hi, b.hi) };
    return z;
}

static inline zmm mulpz4(const zmm& a, const zmm& b) noexcept force_inline;
static inline zmm mulpz4(const zmm& a, const zmm& b) noexcept
{
    const zmm z = { mulpz2(a.lo, b.lo), mulpz2(a.hi, b.hi) };
    return z;
}

template <int N, int mode> static inline zmm scalepz4(const zmm z) noexcept force_inline;
template <int N, int mode> static inline zmm scalepz4(const zmm z) noexcept
{
    constexpr double scale =
        mode == scale_1       ? 1.0           :
        mode == scale_unitary ? 1.0/mysqrt(N) :
        mode == scale_length  ? 1.0/N         : 0.0;
    constexpr ymm sv  = { scale, scale, scale, scale };
    constexpr zmm sv4 = { sv, sv };
    return mode == scale_1 ? z : mulpd4(sv4, z);
}

static inline zmm v8xpz4(const zmm& xy) noexcept force_inline;
static inline zmm v8xpz4(const zmm& xy) noexcept
{
    const zmm z = { v8xpz2(xy.lo), v8xpz2(xy.hi) };
    return z;
}

static inline zmm w8xpz4(const zmm& xy) noexcept force_inline;
static inline zmm w8xpz4(const zmm& xy) noexcept
{
    const zmm z = { w8xpz2(xy.lo), w8xpz2(xy.hi) };
    return z;
}

static inline zmm h1xpz4(const zmm& xy) noexcept force_inline;
static inline zmm h1xpz4(const zmm& xy) noexcept
{
    const zmm z = { h1xpz2(xy.lo), h1xpz2(xy.hi) };
    return z;
}

static inline zmm h3xpz4(const zmm& xy) noexcept force_inline;
static inline zmm h3xpz4(const zmm& xy) noexcept
{
    const zmm z = { h3xpz2(xy.lo), h3xpz2(xy.hi) };
    return z;
}

static inline zmm hfxpz4(const zmm& xy) noexcept force_inline;
static inline zmm hfxpz4(const zmm& xy) noexcept
{
    const zmm z = { hfxpz2(xy.lo), hfxpz2(xy.hi) };
    return z;
}

static inline zmm hdxpz4(const zmm& xy) noexcept force_inline;
static inline zmm hdxpz4(const zmm& xy) noexcept
{
    const zmm z = { hdxpz2(xy.lo), hdxpz2(xy.hi) };
    return z;
}

static inline zmm duppz4(const xmm x) noexcept force_inline;
static inline zmm duppz4(const xmm x) noexcept
{
    const ymm y = duppz2(x);
    const zmm z = { y, y };
    return z;
}

static inline zmm duppz5(const complex_t& z) noexcept force_inline;
static inline zmm duppz5(const complex_t& z) noexcept
{
    return duppz4(getpz(z));
}

} // namespace OTFFT_MISC

#endif // __AVX512F__ && __AVX512DQ__

//=============================================================================
#endif // otfft_misc_h
