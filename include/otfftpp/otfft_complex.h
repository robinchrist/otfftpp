/******************************************************************************
*  OTFFT Complex & Memory Allocator Version 11.4xv
*
*  Copyright (c) 2019 OK Ojisan(Takuya OKAHISA)
*  Released under the MIT license
*  http://opensource.org/licenses/mit-license.php
******************************************************************************/

#ifndef otfft_complex_h
#define otfft_complex_h

#if __GNUC__ >= 3
#define force_inline  __attribute__((const,always_inline))
#define force_inline3 __attribute__((always_inline))
#else
#define force_inline
#define force_inline3
#endif

//=============================================================================
// User Defined Complex Number Class
//=============================================================================

#include <complex>
#include <cstring>
#include <memory>

namespace OTFFT {

struct complex_t
{
    double Re, Im;

    complex_t() noexcept : Re(0), Im(0) {}
    complex_t(double x) noexcept : Re(x), Im(0) {}
    complex_t(double x, double y) noexcept : Re(x), Im(y) {}
    complex_t(const std::complex<double>& z) noexcept : Re(z.real()), Im(z.imag()) {}
    operator std::complex<double>() const { return std::complex<double>(Re, Im); }

    complex_t& operator+=(const complex_t& z) noexcept
    {
        Re += z.Re;
        Im += z.Im;
        return *this;
    }

    complex_t& operator-=(const complex_t& z) noexcept
    {
        Re -= z.Re;
        Im -= z.Im;
        return *this;
    }

    complex_t& operator*=(const double& x) noexcept
    {
        Re *= x;
        Im *= x;
        return *this;
    }

    complex_t& operator/=(const double& x) noexcept
    {
        Re /= x;
        Im /= x;
        return *this;
    }

    complex_t& operator*=(const complex_t& z) noexcept
    {
        const double tmp = Re*z.Re - Im*z.Im;
        Im = Re*z.Im + Im*z.Re;
        Re = tmp;
        return *this;
    }
};

#if 1
typedef double* __restrict const double_vector;
typedef const double* __restrict const const_double_vector;
typedef complex_t* __restrict const complex_vector;
typedef const complex_t* __restrict const const_complex_vector;
#else
typedef double* const double_vector;
typedef const double* const const_double_vector;
typedef complex_t* const complex_vector;
typedef const complex_t* const const_complex_vector;
#endif

static inline double Re(const complex_t& z) noexcept force_inline;
static inline double Re(const complex_t& z) noexcept { return z.Re; }
static inline double Im(const complex_t& z) noexcept force_inline;
static inline double Im(const complex_t& z) noexcept { return z.Im; }

static inline double norm(const complex_t& z) noexcept force_inline;
static inline double norm(const complex_t& z) noexcept
{
    return z.Re*z.Re + z.Im*z.Im;
}
static inline complex_t conj(const complex_t& z) noexcept force_inline;
static inline complex_t conj(const complex_t& z) noexcept
{
    return complex_t(z.Re, -z.Im);
}
static inline double abs(const complex_t& z) noexcept force_inline;
static inline double abs(const complex_t& z) noexcept
{
    return sqrt(norm(z));
}
static inline double arg(const complex_t& z) noexcept force_inline;
static inline double arg(const complex_t& z) noexcept
{
    return atan2(z.Im, z.Re);
}
static inline complex_t proj(const complex_t& z) noexcept force_inline;
static inline complex_t proj(const complex_t& z) noexcept
{
    const double den(norm(z) + double(1.0));
    return complex_t(double(2.0) * z.Re / den, double(2.0) * z.Im / den);
}
static inline complex_t polar(const double rho, const double theta) noexcept force_inline;
static inline complex_t polar(const double rho, const double theta) noexcept
{
    return complex_t(rho * cos(theta), rho * sin(theta));
}
static inline complex_t jx(const complex_t& z) noexcept force_inline;
static inline complex_t jx(const complex_t& z) noexcept
{
    return complex_t(-z.Im, z.Re);
}
static inline complex_t neg(const complex_t& z) noexcept force_inline;
static inline complex_t neg(const complex_t& z) noexcept
{
    return complex_t(-z.Re, -z.Im);
}
static inline complex_t mjx(const complex_t& z) noexcept force_inline;
static inline complex_t mjx(const complex_t& z) noexcept
{
    return complex_t(z.Im, -z.Re);
}
static inline complex_t operator+(const complex_t& a, const complex_t& b) noexcept force_inline;
static inline complex_t operator+(const complex_t& a, const complex_t& b) noexcept
{
    return complex_t(a.Re + b.Re, a.Im + b.Im);
}
static inline complex_t operator-(const complex_t& a, const complex_t& b) noexcept force_inline;
static inline complex_t operator-(const complex_t& a, const complex_t& b) noexcept
{
    return complex_t(a.Re - b.Re, a.Im - b.Im);
}
static inline complex_t operator*(const double& a, const complex_t& b) noexcept force_inline;
static inline complex_t operator*(const double& a, const complex_t& b) noexcept
{
    return complex_t(a*b.Re, a*b.Im);
}
static inline complex_t operator*(const complex_t& a, const complex_t& b) noexcept force_inline;
static inline complex_t operator*(const complex_t& a, const complex_t& b) noexcept
{
    return complex_t(a.Re*b.Re - a.Im*b.Im, a.Re*b.Im + a.Im*b.Re);
}
static inline complex_t operator/(const complex_t& a, const double& b) noexcept force_inline;
static inline complex_t operator/(const complex_t& a, const double& b) noexcept
{
    return complex_t(a.Re/b, a.Im/b);
}
static inline complex_t operator/(const complex_t& a, const complex_t& b) noexcept force_inline;
static inline complex_t operator/(const complex_t& a, const complex_t& b) noexcept
{
    const double b2 = b.Re*b.Re + b.Im*b.Im;
    return (a * conj(b)) / b2;
}

static inline complex_t expj(const double& theta) noexcept force_inline;
static inline complex_t expj(const double& theta) noexcept
{
    //return complex_t(cos(theta), sin(theta));
    return complex_t(std::polar(1.0, theta));
}

} // namespace OTFFT

//=============================================================================
// Aligned Memory Allocator
//=============================================================================

namespace OTFFT {

inline void* generic_aligned_alloc(std::size_t size, std::size_t alignment) noexcept
{
    constexpr size_t N = alignof(void*);

    if (alignment < N) {
        alignment = N;
    }

    std::size_t n = size + alignment - N;
    void* p = std::malloc(sizeof(void*) + n);
    if (p) {
        void* p2 = static_cast<char*>(p) + sizeof(p);
        auto * const ap = static_cast<char*>(std::align(alignment, size, p2, n));
        std::memcpy(ap - sizeof(p), &p, sizeof(p));

        p = ap;
    }
    return p;
}

inline void generic_aligned_free(void* ptr) noexcept
{
    if(ptr) {
        void* rp;
        std::memcpy(&rp, static_cast<char*>(ptr) - sizeof(void*), sizeof(void*));
        std::free(rp);
    }
}


#ifdef __AVX512F__
static inline void* simd_malloc(const size_t n) { return generic_aligned_alloc(n, 64); }
#else
    #ifdef __AVX__
    static inline void* simd_malloc(const size_t n) { return generic_aligned_alloc(n, 32); }
    #else
    static inline void* simd_malloc(const size_t n) { return generic_aligned_alloc(n, 16); }
    #endif
#endif

static inline void simd_free(void* p) { generic_aligned_free(p); }

template <class T> struct simd_array
{
    T* p;

    simd_array() noexcept : p(0) {}
    simd_array(int n) : p((T*) simd_malloc(n*sizeof(T)))
    {
        if (p == 0) throw std::bad_alloc();
    }

    ~simd_array() { if (p) simd_free(p); }

    void setup(int n)
    {
        if (p) simd_free(p);
        p = (T*) simd_malloc(n*sizeof(T));
        if (p == 0) throw std::bad_alloc();
    }

    void destroy() { if (p) simd_free(p); p = 0; }

    T& operator[](int i) noexcept { return p[i]; }
    const T& operator[](int i) const noexcept { return p[i]; }
    T* operator&() const noexcept { return p; }
};

} // namespace OTFFT

//=============================================================================

#endif // otfft_complex_h
