/******************************************************************************
*  OTFFT AVXDIT(Radix-4) of OpenMP Version 11.4xv
*
*  Copyright (c) 2019 OK Ojisan(Takuya OKAHISA)
*  Released under the MIT license
*  http://opensource.org/licenses/mit-license.php
******************************************************************************/

#ifndef otfft_avxdit4omp_h
#define otfft_avxdit4omp_h

#include "otfftpp/otfft_complex.h"
#include "otfftpp/detail/otfft_misc.h"

namespace OTFFT {

namespace OTFFT_AVXDIT4omp { //////////////////////////////////////////////////

    using namespace OTFFT;
    using namespace OTFFT_MISC;

    ///////////////////////////////////////////////////////////////////////////////
    // Forward Butterfly Operation
    ///////////////////////////////////////////////////////////////////////////////

    template <int n, int s> struct fwdcore
    {
        static constexpr int N  = n*s;
        static constexpr int N0 = 0;
        static constexpr int N1 = N/4;
        static constexpr int N2 = N1*2;
        static constexpr int N3 = N1*3;
        static constexpr int Ni = N1/4;
        static constexpr int h  = s/4;

        void operator()(
                complex_vector x, complex_vector y, const_complex_vector W) const noexcept
        {
            #pragma omp for schedule(static)
            for (int i = 0; i < Ni; i++) {
                const int p = i / h;
                const int q = i % h * 4;
                const int sp = s*p;
                const int s4p = 4*sp;
                const zmm w1p = duppz5(W[sp]);
                const zmm w2p = mulpz4(w1p,w1p);
                const zmm w3p = mulpz4(w1p,w2p);
                complex_vector xq_sp  = x + q + sp;
                complex_vector yq_s4p = y + q + s4p;
                const zmm a =             getpz4(yq_s4p+s*0);
                const zmm b = mulpz4(w1p, getpz4(yq_s4p+s*1));
                const zmm c = mulpz4(w2p, getpz4(yq_s4p+s*2));
                const zmm d = mulpz4(w3p, getpz4(yq_s4p+s*3));
                const zmm  apc =       addpz4(a, c);
                const zmm  amc =       subpz4(a, c);
                const zmm  bpd =       addpz4(b, d);
                const zmm jbmd = jxpz4(subpz4(b, d));
                setpz4(xq_sp+N0, addpz4(apc,  bpd));
                setpz4(xq_sp+N1, subpz4(amc, jbmd));
                setpz4(xq_sp+N2, subpz4(apc,  bpd));
                setpz4(xq_sp+N3, addpz4(amc, jbmd));
            }
        }
    };

    template <int N> struct fwdcore<N,1>
    {
        static constexpr int N0 = 0;
        static constexpr int N1 = N/4;
        static constexpr int N2 = N1*2;
        static constexpr int N3 = N1*3;

        void operator()(
                complex_vector x, complex_vector y, const_complex_vector W) const noexcept
        {
            #pragma omp for schedule(static) nowait
            for (int p = 0; p < N1; p += 2) {
                complex_vector x_p  = x + p;
                complex_vector y_4p = y + 4*p;
#if 0
                const ymm w1p = getpz2(W+p);
                const ymm w2p = mulpz2(w1p,w1p);;
                const ymm w3p = mulpz2(w1p,w2p);;
                const ymm a =             getpz3<4>(y_4p+0);
                const ymm b = mulpz2(w1p, getpz3<4>(y_4p+1));
                const ymm c = mulpz2(w2p, getpz3<4>(y_4p+2));
                const ymm d = mulpz2(w3p, getpz3<4>(y_4p+3));
#else
                const ymm w1p = getpz2(W+p);
                const ymm ab = getpz2(y_4p+0);
                const ymm cd = getpz2(y_4p+2);
                const ymm w2p = mulpz2(w1p,w1p);
                const ymm ef = getpz2(y_4p+4);
                const ymm w3p = mulpz2(w1p,w2p);
                const ymm gh = getpz2(y_4p+6);
                const ymm a =             catlo(ab, ef);
                const ymm b = mulpz2(w1p, cathi(ab, ef));
                const ymm c = mulpz2(w2p, catlo(cd, gh));
                const ymm d = mulpz2(w3p, cathi(cd, gh));
#endif
                const ymm  apc =       addpz2(a, c);
                const ymm  amc =       subpz2(a, c);
                const ymm  bpd =       addpz2(b, d);
                const ymm jbmd = jxpz2(subpz2(b, d));
                setpz2(x_p+N0, addpz2(apc,  bpd));
                setpz2(x_p+N1, subpz2(amc, jbmd));
                setpz2(x_p+N2, subpz2(apc,  bpd));
                setpz2(x_p+N3, addpz2(amc, jbmd));
            }
        }
    };

    ///////////////////////////////////////////////////////////////////////////////

    template <int n, int s, bool eo, int mode> struct fwdend;

    //-----------------------------------------------------------------------------

    template <int s, bool eo, int mode> struct fwdend<4,s,eo,mode>
    {
        static constexpr int N = 4*s;

        void operator()(complex_vector x, complex_vector y) const noexcept
        {
            complex_vector z = eo ? y : x;
            #pragma omp for schedule(static)
            for (int q = 0; q < s; q += 2) {
                complex_vector xq = x + q;
                complex_vector zq = z + q;
                const ymm a = scalepz2<N,mode>(getpz2(zq+s*0));
                const ymm b = scalepz2<N,mode>(getpz2(zq+s*1));
                const ymm c = scalepz2<N,mode>(getpz2(zq+s*2));
                const ymm d = scalepz2<N,mode>(getpz2(zq+s*3));
                const ymm  apc =       addpz2(a, c);
                const ymm  amc =       subpz2(a, c);
                const ymm  bpd =       addpz2(b, d);
                const ymm jbmd = jxpz2(subpz2(b, d));
                setpz2(xq+s*0, addpz2(apc,  bpd));
                setpz2(xq+s*1, subpz2(amc, jbmd));
                setpz2(xq+s*2, subpz2(apc,  bpd));
                setpz2(xq+s*3, addpz2(amc, jbmd));
            }
        }
    };

    template <bool eo, int mode> struct fwdend<4,1,eo,mode>
    {
        inline void operator()(complex_vector x, complex_vector y) const noexcept
        {
            #pragma omp single
            {

                complex_vector z = eo ? y : x;
                const xmm a = scalepz<4,mode>(getpz(z[0]));
                const xmm b = scalepz<4,mode>(getpz(z[1]));
                const xmm c = scalepz<4,mode>(getpz(z[2]));
                const xmm d = scalepz<4,mode>(getpz(z[3]));
                const xmm  apc =      addpz(a, c);
                const xmm  amc =      subpz(a, c);
                const xmm  bpd =      addpz(b, d);
                const xmm jbmd = jxpz(subpz(b, d));
                setpz(x[0], addpz(apc,  bpd));
                setpz(x[1], subpz(amc, jbmd));
                setpz(x[2], subpz(apc,  bpd));
                setpz(x[3], addpz(amc, jbmd));
            }
        }
    };

    //-----------------------------------------------------------------------------

    template <int s, bool eo, int mode> struct fwdend<2,s,eo,mode>
    {
        static constexpr int N = 2*s;

        void operator()(complex_vector x, complex_vector y) const noexcept
        {
            complex_vector z = eo ? y : x;
            #pragma omp for schedule(static)
            for (int q = 0; q < s; q += 2) {
                complex_vector xq = x + q;
                complex_vector zq = z + q;
                const ymm a = scalepz2<N,mode>(getpz2(zq+0));
                const ymm b = scalepz2<N,mode>(getpz2(zq+s));
                setpz2(xq+0, addpz2(a, b));
                setpz2(xq+s, subpz2(a, b));
            }
        }
    };

    template <bool eo, int mode> struct fwdend<2,1,eo,mode>
    {
        inline void operator()(complex_vector x, complex_vector y) const noexcept
        {
            #pragma omp single
            {

                complex_vector z = eo ? y : x;
                const xmm a = scalepz<2,mode>(getpz(z[0]));
                const xmm b = scalepz<2,mode>(getpz(z[1]));
                setpz(x[0], addpz(a, b));
                setpz(x[1], subpz(a, b));
            }
        }
    };

    ///////////////////////////////////////////////////////////////////////////////
    // Forward FFT
    ///////////////////////////////////////////////////////////////////////////////

    template <int n, int s, bool eo, int mode> struct fwdfft
    {
        inline void operator()(
                complex_vector x, complex_vector y, const_complex_vector W) const noexcept
        {
            fwdfft<n/4,4*s,!eo,mode>()(y, x, W);
            fwdcore<n,s>()(x, y, W);
        }
    };

    template <int s, bool eo, int mode> struct fwdfft<4,s,eo,mode>
    {
        inline void operator()(
                complex_vector x, complex_vector y, const_complex_vector) const noexcept
        {
            fwdend<4,s,eo,mode>()(x, y);
        }
    };

    template <int s, bool eo, int mode> struct fwdfft<2,s,eo,mode>
    {
        inline void operator()(
                complex_vector x, complex_vector y, const_complex_vector) const noexcept
        {
            fwdend<2,s,eo,mode>()(x, y);
        }
    };

    ///////////////////////////////////////////////////////////////////////////////
    // Inverse Butterfly Operation
    ///////////////////////////////////////////////////////////////////////////////

    template <int n, int s> struct invcore
    {
        static constexpr int N  = n*s;
        static constexpr int N0 = 0;
        static constexpr int N1 = N/4;
        static constexpr int N2 = N1*2;
        static constexpr int N3 = N1*3;
        static constexpr int Ni = N1/4;
        static constexpr int h  = s/4;

        void operator()(
                complex_vector x, complex_vector y, const_complex_vector W) const noexcept
        {
            #pragma omp for schedule(static)
            for (int i = 0; i < Ni; i++) {
                const int p = i / h;
                const int q = i % h * 4;
                const int sp = s*p;
                const int s4p = 4*sp;
                //const zmm w1p = cnjpz4(duppz5(W[sp]));
                const zmm w1p = duppz5(conj(W[sp]));
                const zmm w2p = mulpz4(w1p,w1p);
                const zmm w3p = mulpz4(w1p,w2p);
                complex_vector xq_sp  = x + q + sp;
                complex_vector yq_s4p = y + q + s4p;
                const zmm a =             getpz4(yq_s4p+s*0);
                const zmm b = mulpz4(w1p, getpz4(yq_s4p+s*1));
                const zmm c = mulpz4(w2p, getpz4(yq_s4p+s*2));
                const zmm d = mulpz4(w3p, getpz4(yq_s4p+s*3));
                const zmm  apc =       addpz4(a, c);
                const zmm  amc =       subpz4(a, c);
                const zmm  bpd =       addpz4(b, d);
                const zmm jbmd = jxpz4(subpz4(b, d));
                setpz4(xq_sp+N0, addpz4(apc,  bpd));
                setpz4(xq_sp+N1, addpz4(amc, jbmd));
                setpz4(xq_sp+N2, subpz4(apc,  bpd));
                setpz4(xq_sp+N3, subpz4(amc, jbmd));
            }
        }
    };

    template <int N> struct invcore<N,1>
    {
        static constexpr int N0 = 0;
        static constexpr int N1 = N/4;
        static constexpr int N2 = N1*2;
        static constexpr int N3 = N1*3;

        void operator()(
                complex_vector x, complex_vector y, const_complex_vector W) const noexcept
        {
            #pragma omp for schedule(static) nowait
            for (int p = 0; p < N1; p += 2) {
                complex_vector x_p  = x + p;
                complex_vector y_4p = y + 4*p;
#if 0
                const ymm w1p = cnjpz2(getpz2(W+p));
                const ymm w2p = mulpz2(w1p,w1p);
                const ymm w3p = mulpz2(w1p,w2p);
                const ymm a =             getpz3<4>(y_4p+0);
                const ymm b = mulpz2(w1p, getpz3<4>(y_4p+1));
                const ymm c = mulpz2(w2p, getpz3<4>(y_4p+2));
                const ymm d = mulpz2(w3p, getpz3<4>(y_4p+3));
#else
                const ymm w1p = cnjpz2(getpz2(W+p));
                const ymm ab = getpz2(y_4p+0);
                const ymm cd = getpz2(y_4p+2);
                const ymm w2p = mulpz2(w1p,w1p);
                const ymm ef = getpz2(y_4p+4);
                const ymm w3p = mulpz2(w1p,w2p);
                const ymm gh = getpz2(y_4p+6);
                const ymm a =             catlo(ab, ef);
                const ymm b = mulpz2(w1p, cathi(ab, ef));
                const ymm c = mulpz2(w2p, catlo(cd, gh));
                const ymm d = mulpz2(w3p, cathi(cd, gh));
#endif
                const ymm  apc =       addpz2(a, c);
                const ymm  amc =       subpz2(a, c);
                const ymm  bpd =       addpz2(b, d);
                const ymm jbmd = jxpz2(subpz2(b, d));
                setpz2(x_p+N0, addpz2(apc,  bpd));
                setpz2(x_p+N1, addpz2(amc, jbmd));
                setpz2(x_p+N2, subpz2(apc,  bpd));
                setpz2(x_p+N3, subpz2(amc, jbmd));
            }
        }
    };

    ///////////////////////////////////////////////////////////////////////////////

    template <int n, int s, bool eo, int mode> struct invend;

    //-----------------------------------------------------------------------------

    template <int s, bool eo, int mode> struct invend<4,s,eo,mode>
    {
        static constexpr int N  = 4*s;

        void operator()(complex_vector x, complex_vector y) const noexcept
        {
            complex_vector z = eo ? y : x;
            #pragma omp for schedule(static)
            for (int q = 0; q < s; q += 2) {
                complex_vector xq = x + q;
                complex_vector zq = z + q;
                const ymm a = scalepz2<N,mode>(getpz2(zq+s*0));
                const ymm b = scalepz2<N,mode>(getpz2(zq+s*1));
                const ymm c = scalepz2<N,mode>(getpz2(zq+s*2));
                const ymm d = scalepz2<N,mode>(getpz2(zq+s*3));
                const ymm  apc =       addpz2(a, c);
                const ymm  amc =       subpz2(a, c);
                const ymm  bpd =       addpz2(b, d);
                const ymm jbmd = jxpz2(subpz2(b, d));
                setpz2(xq+s*0, addpz2(apc,  bpd));
                setpz2(xq+s*1, addpz2(amc, jbmd));
                setpz2(xq+s*2, subpz2(apc,  bpd));
                setpz2(xq+s*3, subpz2(amc, jbmd));
            }
        }
    };

    template <bool eo, int mode> struct invend<4,1,eo,mode>
    {
        inline void operator()(complex_vector x, complex_vector y) const noexcept
        {
            #pragma omp single
            {

                complex_vector z = eo ? y : x;
                const xmm a = scalepz<4,mode>(getpz(z[0]));
                const xmm b = scalepz<4,mode>(getpz(z[1]));
                const xmm c = scalepz<4,mode>(getpz(z[2]));
                const xmm d = scalepz<4,mode>(getpz(z[3]));
                const xmm  apc =      addpz(a, c);
                const xmm  amc =      subpz(a, c);
                const xmm  bpd =      addpz(b, d);
                const xmm jbmd = jxpz(subpz(b, d));
                setpz(x[0], addpz(apc,  bpd));
                setpz(x[1], addpz(amc, jbmd));
                setpz(x[2], subpz(apc,  bpd));
                setpz(x[3], subpz(amc, jbmd));
            }
        }
    };

    //-----------------------------------------------------------------------------

    template <int s, bool eo, int mode> struct invend<2,s,eo,mode>
    {
        static constexpr int N  = 2*s;

        void operator()(complex_vector x, complex_vector y) const noexcept
        {
            complex_vector z = eo ? y : x;
            #pragma omp for schedule(static)
            for (int q = 0; q < s; q += 2) {
                complex_vector xq = x + q;
                complex_vector zq = z + q;
                const ymm a = scalepz2<N,mode>(getpz2(zq+0));
                const ymm b = scalepz2<N,mode>(getpz2(zq+s));
                setpz2(xq+0, addpz2(a, b));
                setpz2(xq+s, subpz2(a, b));
            }
        }
    };

    template <bool eo, int mode> struct invend<2,1,eo,mode>
    {
        inline void operator()(complex_vector x, complex_vector y) const noexcept
        {
            #pragma omp single
            {

                complex_vector z = eo ? y : x;
                const xmm a = scalepz<2,mode>(getpz(z[0]));
                const xmm b = scalepz<2,mode>(getpz(z[1]));
                setpz(x[0], addpz(a, b));
                setpz(x[1], subpz(a, b));
            }
        }
    };

    ///////////////////////////////////////////////////////////////////////////////
    // Inverse FFT
    ///////////////////////////////////////////////////////////////////////////////

    template <int n, int s, bool eo, int mode> struct invfft
    {
        inline void operator()(
                complex_vector x, complex_vector y, const_complex_vector W) const noexcept
        {
            invfft<n/4,4*s,!eo,mode>()(y, x, W);
            invcore<n,s>()(x, y, W);
        }
    };

    template <int s, bool eo, int mode> struct invfft<4,s,eo,mode>
    {
        inline void operator()(
                complex_vector x, complex_vector y, const_complex_vector) const noexcept
        {
            invend<4,s,eo,mode>()(x, y);
        }
    };

    template <int s, bool eo, int mode> struct invfft<2,s,eo,mode>
    {
        inline void operator()(
                complex_vector x, complex_vector y, const_complex_vector) const noexcept
        {
            invend<2,s,eo,mode>()(x, y);
        }
    };

    ///////////////////////////////////////////////////////////////////////////////
    // Power of 2 FFT Routine
    ///////////////////////////////////////////////////////////////////////////////

    inline void fwd(const int log_N,
                    complex_vector x, complex_vector y, const_complex_vector W) noexcept
    {
        constexpr int mode = scale_length;
#pragma omp parallel firstprivate(x,y,W)
        switch (log_N) {
        case  0: break;
        case  1: fwdfft<(1<< 1),1,0,mode>()(x, y, W); break;
        case  2: fwdfft<(1<< 2),1,0,mode>()(x, y, W); break;
        case  3: fwdfft<(1<< 3),1,0,mode>()(x, y, W); break;
        case  4: fwdfft<(1<< 4),1,0,mode>()(x, y, W); break;
        case  5: fwdfft<(1<< 5),1,0,mode>()(x, y, W); break;
        case  6: fwdfft<(1<< 6),1,0,mode>()(x, y, W); break;
        case  7: fwdfft<(1<< 7),1,0,mode>()(x, y, W); break;
        case  8: fwdfft<(1<< 8),1,0,mode>()(x, y, W); break;
        case  9: fwdfft<(1<< 9),1,0,mode>()(x, y, W); break;
        case 10: fwdfft<(1<<10),1,0,mode>()(x, y, W); break;
        case 11: fwdfft<(1<<11),1,0,mode>()(x, y, W); break;
        case 12: fwdfft<(1<<12),1,0,mode>()(x, y, W); break;
        case 13: fwdfft<(1<<13),1,0,mode>()(x, y, W); break;
        case 14: fwdfft<(1<<14),1,0,mode>()(x, y, W); break;
        case 15: fwdfft<(1<<15),1,0,mode>()(x, y, W); break;
        case 16: fwdfft<(1<<16),1,0,mode>()(x, y, W); break;
        case 17: fwdfft<(1<<17),1,0,mode>()(x, y, W); break;
        case 18: fwdfft<(1<<18),1,0,mode>()(x, y, W); break;
        case 19: fwdfft<(1<<19),1,0,mode>()(x, y, W); break;
        case 20: fwdfft<(1<<20),1,0,mode>()(x, y, W); break;
        case 21: fwdfft<(1<<21),1,0,mode>()(x, y, W); break;
        case 22: fwdfft<(1<<22),1,0,mode>()(x, y, W); break;
        case 23: fwdfft<(1<<23),1,0,mode>()(x, y, W); break;
        case 24: fwdfft<(1<<24),1,0,mode>()(x, y, W); break;
        }
    }

    inline void fwd0(const int log_N,
                     complex_vector x, complex_vector y, const_complex_vector W) noexcept
    {
        constexpr int mode = scale_1;
#pragma omp parallel firstprivate(x,y,W)
        switch (log_N) {
        case  0: break;
        case  1: fwdfft<(1<< 1),1,0,mode>()(x, y, W); break;
        case  2: fwdfft<(1<< 2),1,0,mode>()(x, y, W); break;
        case  3: fwdfft<(1<< 3),1,0,mode>()(x, y, W); break;
        case  4: fwdfft<(1<< 4),1,0,mode>()(x, y, W); break;
        case  5: fwdfft<(1<< 5),1,0,mode>()(x, y, W); break;
        case  6: fwdfft<(1<< 6),1,0,mode>()(x, y, W); break;
        case  7: fwdfft<(1<< 7),1,0,mode>()(x, y, W); break;
        case  8: fwdfft<(1<< 8),1,0,mode>()(x, y, W); break;
        case  9: fwdfft<(1<< 9),1,0,mode>()(x, y, W); break;
        case 10: fwdfft<(1<<10),1,0,mode>()(x, y, W); break;
        case 11: fwdfft<(1<<11),1,0,mode>()(x, y, W); break;
        case 12: fwdfft<(1<<12),1,0,mode>()(x, y, W); break;
        case 13: fwdfft<(1<<13),1,0,mode>()(x, y, W); break;
        case 14: fwdfft<(1<<14),1,0,mode>()(x, y, W); break;
        case 15: fwdfft<(1<<15),1,0,mode>()(x, y, W); break;
        case 16: fwdfft<(1<<16),1,0,mode>()(x, y, W); break;
        case 17: fwdfft<(1<<17),1,0,mode>()(x, y, W); break;
        case 18: fwdfft<(1<<18),1,0,mode>()(x, y, W); break;
        case 19: fwdfft<(1<<19),1,0,mode>()(x, y, W); break;
        case 20: fwdfft<(1<<20),1,0,mode>()(x, y, W); break;
        case 21: fwdfft<(1<<21),1,0,mode>()(x, y, W); break;
        case 22: fwdfft<(1<<22),1,0,mode>()(x, y, W); break;
        case 23: fwdfft<(1<<23),1,0,mode>()(x, y, W); break;
        case 24: fwdfft<(1<<24),1,0,mode>()(x, y, W); break;
        }
    }

    inline void fwdu(const int log_N,
                     complex_vector x, complex_vector y, const_complex_vector W) noexcept
    {
        constexpr int mode = scale_unitary;
#pragma omp parallel firstprivate(x,y,W)
        switch (log_N) {
        case  0: break;
        case  1: fwdfft<(1<< 1),1,0,mode>()(x, y, W); break;
        case  2: fwdfft<(1<< 2),1,0,mode>()(x, y, W); break;
        case  3: fwdfft<(1<< 3),1,0,mode>()(x, y, W); break;
        case  4: fwdfft<(1<< 4),1,0,mode>()(x, y, W); break;
        case  5: fwdfft<(1<< 5),1,0,mode>()(x, y, W); break;
        case  6: fwdfft<(1<< 6),1,0,mode>()(x, y, W); break;
        case  7: fwdfft<(1<< 7),1,0,mode>()(x, y, W); break;
        case  8: fwdfft<(1<< 8),1,0,mode>()(x, y, W); break;
        case  9: fwdfft<(1<< 9),1,0,mode>()(x, y, W); break;
        case 10: fwdfft<(1<<10),1,0,mode>()(x, y, W); break;
        case 11: fwdfft<(1<<11),1,0,mode>()(x, y, W); break;
        case 12: fwdfft<(1<<12),1,0,mode>()(x, y, W); break;
        case 13: fwdfft<(1<<13),1,0,mode>()(x, y, W); break;
        case 14: fwdfft<(1<<14),1,0,mode>()(x, y, W); break;
        case 15: fwdfft<(1<<15),1,0,mode>()(x, y, W); break;
        case 16: fwdfft<(1<<16),1,0,mode>()(x, y, W); break;
        case 17: fwdfft<(1<<17),1,0,mode>()(x, y, W); break;
        case 18: fwdfft<(1<<18),1,0,mode>()(x, y, W); break;
        case 19: fwdfft<(1<<19),1,0,mode>()(x, y, W); break;
        case 20: fwdfft<(1<<20),1,0,mode>()(x, y, W); break;
        case 21: fwdfft<(1<<21),1,0,mode>()(x, y, W); break;
        case 22: fwdfft<(1<<22),1,0,mode>()(x, y, W); break;
        case 23: fwdfft<(1<<23),1,0,mode>()(x, y, W); break;
        case 24: fwdfft<(1<<24),1,0,mode>()(x, y, W); break;
        }
    }

    inline void fwdn(const int log_N,
                     complex_vector x, complex_vector y, const_complex_vector W) noexcept
    {
        fwd(log_N, x, y, W);
    }

    ///////////////////////////////////////////////////////////////////////////////

    inline void inv(const int log_N,
                    complex_vector x, complex_vector y, const_complex_vector W) noexcept
    {
        constexpr int mode = scale_1;
#pragma omp parallel firstprivate(x,y,W)
        switch (log_N) {
        case  0: break;
        case  1: invfft<(1<< 1),1,0,mode>()(x, y, W); break;
        case  2: invfft<(1<< 2),1,0,mode>()(x, y, W); break;
        case  3: invfft<(1<< 3),1,0,mode>()(x, y, W); break;
        case  4: invfft<(1<< 4),1,0,mode>()(x, y, W); break;
        case  5: invfft<(1<< 5),1,0,mode>()(x, y, W); break;
        case  6: invfft<(1<< 6),1,0,mode>()(x, y, W); break;
        case  7: invfft<(1<< 7),1,0,mode>()(x, y, W); break;
        case  8: invfft<(1<< 8),1,0,mode>()(x, y, W); break;
        case  9: invfft<(1<< 9),1,0,mode>()(x, y, W); break;
        case 10: invfft<(1<<10),1,0,mode>()(x, y, W); break;
        case 11: invfft<(1<<11),1,0,mode>()(x, y, W); break;
        case 12: invfft<(1<<12),1,0,mode>()(x, y, W); break;
        case 13: invfft<(1<<13),1,0,mode>()(x, y, W); break;
        case 14: invfft<(1<<14),1,0,mode>()(x, y, W); break;
        case 15: invfft<(1<<15),1,0,mode>()(x, y, W); break;
        case 16: invfft<(1<<16),1,0,mode>()(x, y, W); break;
        case 17: invfft<(1<<17),1,0,mode>()(x, y, W); break;
        case 18: invfft<(1<<18),1,0,mode>()(x, y, W); break;
        case 19: invfft<(1<<19),1,0,mode>()(x, y, W); break;
        case 20: invfft<(1<<20),1,0,mode>()(x, y, W); break;
        case 21: invfft<(1<<21),1,0,mode>()(x, y, W); break;
        case 22: invfft<(1<<22),1,0,mode>()(x, y, W); break;
        case 23: invfft<(1<<23),1,0,mode>()(x, y, W); break;
        case 24: invfft<(1<<24),1,0,mode>()(x, y, W); break;
        }
    }

    inline void inv0(const int log_N,
                     complex_vector x, complex_vector y, const_complex_vector W) noexcept
    {
        inv(log_N, x, y, W);
    }

    inline void invu(const int log_N,
                     complex_vector x, complex_vector y, const_complex_vector W) noexcept
    {
        constexpr int mode = scale_unitary;
#pragma omp parallel firstprivate(x,y,W)
        switch (log_N) {
        case  0: break;
        case  1: invfft<(1<< 1),1,0,mode>()(x, y, W); break;
        case  2: invfft<(1<< 2),1,0,mode>()(x, y, W); break;
        case  3: invfft<(1<< 3),1,0,mode>()(x, y, W); break;
        case  4: invfft<(1<< 4),1,0,mode>()(x, y, W); break;
        case  5: invfft<(1<< 5),1,0,mode>()(x, y, W); break;
        case  6: invfft<(1<< 6),1,0,mode>()(x, y, W); break;
        case  7: invfft<(1<< 7),1,0,mode>()(x, y, W); break;
        case  8: invfft<(1<< 8),1,0,mode>()(x, y, W); break;
        case  9: invfft<(1<< 9),1,0,mode>()(x, y, W); break;
        case 10: invfft<(1<<10),1,0,mode>()(x, y, W); break;
        case 11: invfft<(1<<11),1,0,mode>()(x, y, W); break;
        case 12: invfft<(1<<12),1,0,mode>()(x, y, W); break;
        case 13: invfft<(1<<13),1,0,mode>()(x, y, W); break;
        case 14: invfft<(1<<14),1,0,mode>()(x, y, W); break;
        case 15: invfft<(1<<15),1,0,mode>()(x, y, W); break;
        case 16: invfft<(1<<16),1,0,mode>()(x, y, W); break;
        case 17: invfft<(1<<17),1,0,mode>()(x, y, W); break;
        case 18: invfft<(1<<18),1,0,mode>()(x, y, W); break;
        case 19: invfft<(1<<19),1,0,mode>()(x, y, W); break;
        case 20: invfft<(1<<20),1,0,mode>()(x, y, W); break;
        case 21: invfft<(1<<21),1,0,mode>()(x, y, W); break;
        case 22: invfft<(1<<22),1,0,mode>()(x, y, W); break;
        case 23: invfft<(1<<23),1,0,mode>()(x, y, W); break;
        case 24: invfft<(1<<24),1,0,mode>()(x, y, W); break;
        }
    }

    inline void invn(const int log_N,
                     complex_vector x, complex_vector y, const_complex_vector W) noexcept
    {
        constexpr int mode = scale_length;
#pragma omp parallel firstprivate(x,y,W)
        switch (log_N) {
        case  0: break;
        case  1: invfft<(1<< 1),1,0,mode>()(x, y, W); break;
        case  2: invfft<(1<< 2),1,0,mode>()(x, y, W); break;
        case  3: invfft<(1<< 3),1,0,mode>()(x, y, W); break;
        case  4: invfft<(1<< 4),1,0,mode>()(x, y, W); break;
        case  5: invfft<(1<< 5),1,0,mode>()(x, y, W); break;
        case  6: invfft<(1<< 6),1,0,mode>()(x, y, W); break;
        case  7: invfft<(1<< 7),1,0,mode>()(x, y, W); break;
        case  8: invfft<(1<< 8),1,0,mode>()(x, y, W); break;
        case  9: invfft<(1<< 9),1,0,mode>()(x, y, W); break;
        case 10: invfft<(1<<10),1,0,mode>()(x, y, W); break;
        case 11: invfft<(1<<11),1,0,mode>()(x, y, W); break;
        case 12: invfft<(1<<12),1,0,mode>()(x, y, W); break;
        case 13: invfft<(1<<13),1,0,mode>()(x, y, W); break;
        case 14: invfft<(1<<14),1,0,mode>()(x, y, W); break;
        case 15: invfft<(1<<15),1,0,mode>()(x, y, W); break;
        case 16: invfft<(1<<16),1,0,mode>()(x, y, W); break;
        case 17: invfft<(1<<17),1,0,mode>()(x, y, W); break;
        case 18: invfft<(1<<18),1,0,mode>()(x, y, W); break;
        case 19: invfft<(1<<19),1,0,mode>()(x, y, W); break;
        case 20: invfft<(1<<20),1,0,mode>()(x, y, W); break;
        case 21: invfft<(1<<21),1,0,mode>()(x, y, W); break;
        case 22: invfft<(1<<22),1,0,mode>()(x, y, W); break;
        case 23: invfft<(1<<23),1,0,mode>()(x, y, W); break;
        case 24: invfft<(1<<24),1,0,mode>()(x, y, W); break;
        }
    }

} /////////////////////////////////////////////////////////////////////////////

}

#endif // otfft_avxdit4omp_h
