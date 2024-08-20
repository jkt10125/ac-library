#include <bits/stdc++.h>

#ifndef ATCODER_MATRIX_HPP
#define ATCODER_MATRIX_HPP 1

#include <vector>
#include <array>
#include <cassert>

namespace atcoder {

template <class T> struct matrix {
  public:
    matrix() : rows_count(0), columns_count(0) {}
    matrix(unsigned int n, unsigned int m) : d(n, std::vector<T>(m)), rows_count(n), columns_count(m) { }

    matrix(std::vector<std::vector<T>> d) : d(std::move(d)), rows_count(this->d.size()), columns_count(this->d[0].size()) { }

    std::vector<T>& operator[](unsigned int i) { return d[i]; }
    const std::vector<T>& operator[](unsigned int i) const { return d[i]; }

    static matrix identity(unsigned int n) {
        matrix mat(n, n);
        for (int i = 0; i < int(n); i++) mat[i][i] = 1;
        return mat;
    }

    matrix transpose() const {
        matrix r(columns_count, rows_count);
        for (unsigned int i = 0; i < rows_count; i++) {
            for (unsigned int j = 0; j < columns_count; j++) {
                r[j][i] = d[i][j];
            }
        }
        return r;
    }

    matrix operator-() const {
        matrix r(rows_count, columns_count);
        for (unsigned int i = 0; i < rows_count; i++) {
            for (unsigned int j = 0; j < columns_count; j++) {
                r[i][j] = -d[i][j];
            }
        }
        return r;
    }

    matrix& operator+=(const matrix& r) {
        for (unsigned int i = 0; i < rows_count; i++) {
            for (unsigned int j = 0; j < columns_count; j++) {
                d[i][j] += r[i][j];
            }
        }
        return *this;
    }

    matrix operator+(const matrix& r) const {
        return matrix(*this) += r;
    }

    matrix& operator-=(const matrix& r) {
        return *this += -r;
    }

    matrix operator-(const matrix& r) const {
        return matrix(*this) -= r;
    }

    matrix operator*(matrix r) const {
        assert(columns_count == r.rows_count);
        r = std::move(r.transpose());
        matrix res(rows_count, r.columns_count);
        for (unsigned int i = 0; i < rows_count; i++) {
            for (unsigned int j = 0; j < r.columns_count; j++) {
                T sum = 0;
                for (unsigned int k = 0; k < columns_count; k++) {
                    sum += d[i][k] * r[j][k];
                }
                res[i][j] = sum;
            }
        }
        return res;
    }

    matrix& operator*=(const matrix& r) {
        return *this = *this * r;
    }

    matrix& operator*=(const T& r) {
        for (unsigned int i = 0; i < rows_count; i++) {
            for (unsigned int j = 0; j < columns_count; j++) {
                d[i][j] *= r;
            }
        }
        return *this;
    }

    matrix operator*(const T& r) const {
        return matrix(*this) *= r;
    }

    matrix& operator/=(const T& r) {
        for (unsigned int i = 0; i < rows_count; i++) {
            for (unsigned int j = 0; j < columns_count; j++) {
                d[i][j] /= r;
            }
        }
        return *this;
    }

    matrix operator/(const T& r) const {
        return matrix(*this) /= r;
    }

    matrix pow(long long k) const {
        matrix res = identity(rows_count);
        for (matrix a = (*this); k; k >>= 1, a *= a) {
            if (k & 1) res *= a;
        }
        return res;
    }

    T determinant() const {
        matrix a = *this;
        T res = T(1);
        for (unsigned int i = 0; i < rows_count; ++i) {
            unsigned int idx = -1;
            for (unsigned int j = i; j < rows_count; ++j) {
                if (a[j][i] != T(0)) {
                    idx = j;
                    break;
                }
            }
            if (idx == -1) {
                return T(0);
            }
            if (idx != i) {
                res = -res;
                for (unsigned int j = i; j < rows_count; ++j) {
                    std::swap(a[i][j], a[idx][j]);
                }
            }
            res *= a[i][i];
            T inv = T(1) / a[i][i];
            for (unsigned int j = i; j < rows_count; ++j) {
                a[i][j] *= inv;
            }
            for (unsigned int j = i + 1; j < rows_count; ++j) {
                T mul = a[j][i];
                for (unsigned int k = i; k < rows_count; ++k) {
                    a[j][k] -= a[i][k] * mul;
                }
            }
        }
        return res;
    }

    matrix inverse() const {
        matrix a = *this;
        std::vector<std::array<unsigned int, 2>> swaps;
        for (unsigned int i = 0; i < rows_count; ++i) {
            unsigned int idx = -1;
            for (unsigned int j = i; j < rows_count; ++j) {
                if (a[j][i] != T(0)) {
                    idx = j;
                    break;
                }
            }
            if (idx == -1) {
                return {};
            }
            if (idx != i) {
                swaps.push_back({idx, i});
                a[idx].swap(a[i]);
            }
            a[i][i] = T(1) / a[i][i];
            for (unsigned int j = 0; j < rows_count; ++j) {
                if (j == i) continue;
                a[i][j] *= a[i][i];
            }
            for (unsigned int j = 0; j < rows_count; ++j) {
                if (j == i) continue;
                for (unsigned int k = 0; k < rows_count; ++k) {
                    if (k == i) continue;
                    a[j][k] -= a[j][i] * a[i][k];
                }
                a[j][i] *= -a[i][i];
            }
        }
        for (unsigned int i = unsigned(swaps.size()); i--;) {
            for (unsigned int j = 0; j < rows_count; ++j) {
                std::swap(a[j][swaps[i][0]], a[j][swaps[i][1]]);
            }
        }
        return a;
    }

//   private:
    std::vector<std::vector<T>> d;
    unsigned int rows_count, columns_count;
};

}  // namespace atcoder

#endif  // ATCODER_MATRIX_HPP

#ifndef ATCODER_COMBINATORICS_HPP
#define ATCODER_COMBINATORICS_HPP 1

#ifndef ATCODER_MODINT_HPP
#define ATCODER_MODINT_HPP 1

#ifndef ATCODER_INTERNAL_MATH_HPP
#define ATCODER_INTERNAL_MATH_HPP 1

#include <utility>

#ifdef _MSC_VER
#include <intrin.h>
#endif

namespace atcoder {

namespace internal {

// @param m `1 <= m`
// @return x mod m
constexpr long long safe_mod(long long x, long long m) {
    x %= m;
    if (x < 0) x += m;
    return x;
}

// Fast modular multiplication by barrett reduction
// Reference: https://en.wikipedia.org/wiki/Barrett_reduction
// NOTE: reconsider after Ice Lake
struct barrett {
    unsigned int _m;
    unsigned long long im;

    // @param m `1 <= m`
    explicit barrett(unsigned int m) : _m(m), im((unsigned long long)(-1) / m + 1) {}

    // @return m
    unsigned int umod() const { return _m; }

    // @param a `0 <= a < m`
    // @param b `0 <= b < m`
    // @return `a * b % m`
    unsigned int mul(unsigned int a, unsigned int b) const {
        // [1] m = 1
        // a = b = im = 0, so okay

        // [2] m >= 2
        // im = ceil(2^64 / m)
        // -> im * m = 2^64 + r (0 <= r < m)
        // let z = a*b = c*m + d (0 <= c, d < m)
        // a*b * im = (c*m + d) * im = c*(im*m) + d*im = c*2^64 + c*r + d*im
        // c*r + d*im < m * m + m * im < m * m + 2^64 + m <= 2^64 + m * (m + 1) < 2^64 * 2
        // ((ab * im) >> 64) == c or c + 1
        unsigned long long z = a;
        z *= b;
#ifdef _MSC_VER
        unsigned long long x;
        _umul128(z, im, &x);
#else
        unsigned long long x =
            (unsigned long long)(((unsigned __int128)(z)*im) >> 64);
#endif
        unsigned long long y = x * _m;
        return (unsigned int)(z - y + (z < y ? _m : 0));
    }
};

// @param n `0 <= n`
// @param m `1 <= m`
// @return `(x ** n) % m`
constexpr long long pow_mod_constexpr(long long x, long long n, int m) {
    if (m == 1) return 0;
    unsigned int _m = (unsigned int)(m);
    unsigned long long r = 1;
    unsigned long long y = safe_mod(x, m);
    while (n) {
        if (n & 1) r = (r * y) % _m;
        y = (y * y) % _m;
        n >>= 1;
    }
    return r;
}

// Reference:
// M. Forisek and J. Jancina,
// Fast Primality Testing for Integers That Fit into a Machine Word
// @param n `0 <= n`
constexpr bool is_prime_constexpr(int n) {
    if (n <= 1) return false;
    if (n == 2 || n == 7 || n == 61) return true;
    if (n % 2 == 0) return false;
    long long d = n - 1;
    while (d % 2 == 0) d /= 2;
    constexpr long long bases[3] = {2, 7, 61};
    for (long long a : bases) {
        long long t = d;
        long long y = pow_mod_constexpr(a, t, n);
        while (t != n - 1 && y != 1 && y != n - 1) {
            y = y * y % n;
            t <<= 1;
        }
        if (y != n - 1 && t % 2 == 0) {
            return false;
        }
    }
    return true;
}
template <int n> constexpr bool is_prime = is_prime_constexpr(n);

// @param b `1 <= b`
// @return pair(g, x) s.t. g = gcd(a, b), xa = g (mod b), 0 <= x < b/g
constexpr std::pair<long long, long long> inv_gcd(long long a, long long b) {
    a = safe_mod(a, b);
    if (a == 0) return {b, 0};

    // Contracts:
    // [1] s - m0 * a = 0 (mod b)
    // [2] t - m1 * a = 0 (mod b)
    // [3] s * |m1| + t * |m0| <= b
    long long s = b, t = a;
    long long m0 = 0, m1 = 1;

    while (t) {
        long long u = s / t;
        s -= t * u;
        m0 -= m1 * u;  // |m1 * u| <= |m1| * s <= b

        // [3]:
        // (s - t * u) * |m1| + t * |m0 - m1 * u|
        // <= s * |m1| - t * u * |m1| + t * (|m0| + |m1| * u)
        // = s * |m1| + t * |m0| <= b

        auto tmp = s;
        s = t;
        t = tmp;
        tmp = m0;
        m0 = m1;
        m1 = tmp;
    }
    // by [3]: |m0| <= b/g
    // by g != b: |m0| < b/g
    if (m0 < 0) m0 += b / s;
    return {s, m0};
}

// Compile time primitive root
// @param m must be prime
// @return primitive root (and minimum in now)
constexpr int primitive_root_constexpr(int m) {
    if (m == 2) return 1;
    if (m == 167772161) return 3;
    if (m == 469762049) return 3;
    if (m == 754974721) return 11;
    if (m == 998244353) return 3;
    int divs[20] = {};
    divs[0] = 2;
    int cnt = 1;
    int x = (m - 1) / 2;
    while (x % 2 == 0) x /= 2;
    for (int i = 3; (long long)(i)*i <= x; i += 2) {
        if (x % i == 0) {
            divs[cnt++] = i;
            while (x % i == 0) {
                x /= i;
            }
        }
    }
    if (x > 1) {
        divs[cnt++] = x;
    }
    for (int g = 2;; g++) {
        bool ok = true;
        for (int i = 0; i < cnt; i++) {
            if (pow_mod_constexpr(g, (m - 1) / divs[i], m) == 1) {
                ok = false;
                break;
            }
        }
        if (ok) return g;
    }
}
template <int m> constexpr int primitive_root = primitive_root_constexpr(m);

// @param n `n < 2^32`
// @param m `1 <= m < 2^32`
// @return sum_{i=0}^{n-1} floor((ai + b) / m) (mod 2^64)
unsigned long long floor_sum_unsigned(unsigned long long n,
                                      unsigned long long m,
                                      unsigned long long a,
                                      unsigned long long b) {
    unsigned long long ans = 0;
    while (true) {
        if (a >= m) {
            ans += n * (n - 1) / 2 * (a / m);
            a %= m;
        }
        if (b >= m) {
            ans += n * (b / m);
            b %= m;
        }

        unsigned long long y_max = a * n + b;
        if (y_max < m) break;
        // y_max < m * (n + 1)
        // floor(y_max / m) <= n
        n = (unsigned long long)(y_max / m);
        b = (unsigned long long)(y_max % m);
        std::swap(m, a);
    }
    return ans;
}

}  // namespace internal

}  // namespace atcoder

#endif  // ATCODER_INTERNAL_MATH_HPP

#ifndef ATCODER_INTERNAL_TYPE_TRAITS_HPP
#define ATCODER_INTERNAL_TYPE_TRAITS_HPP 1

#include <cassert>
#include <numeric>
#include <type_traits>

namespace atcoder {

namespace internal {

#ifndef _MSC_VER
template <class T>
using is_signed_int128 =
    typename std::conditional<std::is_same<T, __int128_t>::value ||
                                  std::is_same<T, __int128>::value,
                              std::true_type,
                              std::false_type>::type;

template <class T>
using is_unsigned_int128 =
    typename std::conditional<std::is_same<T, __uint128_t>::value ||
                                  std::is_same<T, unsigned __int128>::value,
                              std::true_type,
                              std::false_type>::type;

template <class T>
using make_unsigned_int128 =
    typename std::conditional<std::is_same<T, __int128_t>::value,
                              __uint128_t,
                              unsigned __int128>;

template <class T>
using is_integral = typename std::conditional<std::is_integral<T>::value ||
                                                  is_signed_int128<T>::value ||
                                                  is_unsigned_int128<T>::value,
                                              std::true_type,
                                              std::false_type>::type;

template <class T>
using is_signed_int = typename std::conditional<(is_integral<T>::value &&
                                                 std::is_signed<T>::value) ||
                                                    is_signed_int128<T>::value,
                                                std::true_type,
                                                std::false_type>::type;

template <class T>
using is_unsigned_int =
    typename std::conditional<(is_integral<T>::value &&
                               std::is_unsigned<T>::value) ||
                                  is_unsigned_int128<T>::value,
                              std::true_type,
                              std::false_type>::type;

template <class T>
using to_unsigned = typename std::conditional<
    is_signed_int128<T>::value,
    make_unsigned_int128<T>,
    typename std::conditional<std::is_signed<T>::value,
                              std::make_unsigned<T>,
                              std::common_type<T>>::type>::type;

#else

template <class T> using is_integral = typename std::is_integral<T>;

template <class T>
using is_signed_int =
    typename std::conditional<is_integral<T>::value && std::is_signed<T>::value,
                              std::true_type,
                              std::false_type>::type;

template <class T>
using is_unsigned_int =
    typename std::conditional<is_integral<T>::value &&
                                  std::is_unsigned<T>::value,
                              std::true_type,
                              std::false_type>::type;

template <class T>
using to_unsigned = typename std::conditional<is_signed_int<T>::value,
                                              std::make_unsigned<T>,
                                              std::common_type<T>>::type;

#endif

template <class T>
using is_signed_int_t = std::enable_if_t<is_signed_int<T>::value>;

template <class T>
using is_unsigned_int_t = std::enable_if_t<is_unsigned_int<T>::value>;

template <class T> using to_unsigned_t = typename to_unsigned<T>::type;

}  // namespace internal

}  // namespace atcoder

#endif  // ATCODER_INTERNAL_TYPE_TRAITS_HPP


#include <cassert>
#include <numeric>
#include <type_traits>

#ifdef _MSC_VER
#include <intrin.h>
#endif

namespace atcoder {

namespace internal {

struct modint_base {};
struct static_modint_base : modint_base {};

template <class T> using is_modint = std::is_base_of<modint_base, T>;
template <class T> using is_modint_t = std::enable_if_t<is_modint<T>::value>;

}  // namespace internal

template <int m, std::enable_if_t<(1 <= m)>* = nullptr>
struct static_modint : internal::static_modint_base {
    using mint = static_modint;

  public:
    static constexpr int mod() { return m; }
    static mint raw(int v) {
        mint x;
        x._v = v;
        return x;
    }

    static_modint() : _v(0) {}
    template <class T, internal::is_signed_int_t<T>* = nullptr>
    static_modint(T v) {
        long long x = (long long)(v % (long long)(umod()));
        if (x < 0) x += umod();
        _v = (unsigned int)(x);
    }
    template <class T, internal::is_unsigned_int_t<T>* = nullptr>
    static_modint(T v) {
        _v = (unsigned int)(v % umod());
    }

    unsigned int val() const { return _v; }

    mint& operator++() {
        _v++;
        if (_v == umod()) _v = 0;
        return *this;
    }
    mint& operator--() {
        if (_v == 0) _v = umod();
        _v--;
        return *this;
    }
    mint operator++(int) {
        mint result = *this;
        ++*this;
        return result;
    }
    mint operator--(int) {
        mint result = *this;
        --*this;
        return result;
    }

    mint& operator+=(const mint& rhs) {
        _v += rhs._v;
        if (_v >= umod()) _v -= umod();
        return *this;
    }
    mint& operator-=(const mint& rhs) {
        _v -= rhs._v;
        if (_v >= umod()) _v += umod();
        return *this;
    }
    mint& operator*=(const mint& rhs) {
        unsigned long long z = _v;
        z *= rhs._v;
        _v = (unsigned int)(z % umod());
        return *this;
    }
    mint& operator/=(const mint& rhs) { return *this = *this * rhs.inv(); }

    mint operator+() const { return *this; }
    mint operator-() const { return mint() - *this; }

    mint pow(long long n) const {
        assert(0 <= n);
        mint x = *this, r = 1;
        while (n) {
            if (n & 1) r *= x;
            x *= x;
            n >>= 1;
        }
        return r;
    }
    mint inv() const {
        if (prime) {
            assert(_v);
            return pow(umod() - 2);
        } else {
            auto eg = internal::inv_gcd(_v, m);
            assert(eg.first == 1);
            return eg.second;
        }
    }

    friend mint operator+(const mint& lhs, const mint& rhs) {
        return mint(lhs) += rhs;
    }
    friend mint operator-(const mint& lhs, const mint& rhs) {
        return mint(lhs) -= rhs;
    }
    friend mint operator*(const mint& lhs, const mint& rhs) {
        return mint(lhs) *= rhs;
    }
    friend mint operator/(const mint& lhs, const mint& rhs) {
        return mint(lhs) /= rhs;
    }
    friend bool operator==(const mint& lhs, const mint& rhs) {
        return lhs._v == rhs._v;
    }
    friend bool operator!=(const mint& lhs, const mint& rhs) {
        return lhs._v != rhs._v;
    }

  private:
    unsigned int _v;
    static constexpr unsigned int umod() { return m; }
    static constexpr bool prime = internal::is_prime<m>;
};

template <int id> struct dynamic_modint : internal::modint_base {
    using mint = dynamic_modint;

  public:
    static int mod() { return (int)(bt.umod()); }
    static void set_mod(int m) {
        assert(1 <= m);
        bt = internal::barrett(m);
    }
    static mint raw(int v) {
        mint x;
        x._v = v;
        return x;
    }

    dynamic_modint() : _v(0) {}
    template <class T, internal::is_signed_int_t<T>* = nullptr>
    dynamic_modint(T v) {
        long long x = (long long)(v % (long long)(mod()));
        if (x < 0) x += mod();
        _v = (unsigned int)(x);
    }
    template <class T, internal::is_unsigned_int_t<T>* = nullptr>
    dynamic_modint(T v) {
        _v = (unsigned int)(v % mod());
    }

    unsigned int val() const { return _v; }

    mint& operator++() {
        _v++;
        if (_v == umod()) _v = 0;
        return *this;
    }
    mint& operator--() {
        if (_v == 0) _v = umod();
        _v--;
        return *this;
    }
    mint operator++(int) {
        mint result = *this;
        ++*this;
        return result;
    }
    mint operator--(int) {
        mint result = *this;
        --*this;
        return result;
    }

    mint& operator+=(const mint& rhs) {
        _v += rhs._v;
        if (_v >= umod()) _v -= umod();
        return *this;
    }
    mint& operator-=(const mint& rhs) {
        _v += mod() - rhs._v;
        if (_v >= umod()) _v -= umod();
        return *this;
    }
    mint& operator*=(const mint& rhs) {
        _v = bt.mul(_v, rhs._v);
        return *this;
    }
    mint& operator/=(const mint& rhs) { return *this = *this * rhs.inv(); }

    mint operator+() const { return *this; }
    mint operator-() const { return mint() - *this; }

    mint pow(long long n) const {
        assert(0 <= n);
        mint x = *this, r = 1;
        while (n) {
            if (n & 1) r *= x;
            x *= x;
            n >>= 1;
        }
        return r;
    }
    mint inv() const {
        auto eg = internal::inv_gcd(_v, mod());
        assert(eg.first == 1);
        return eg.second;
    }

    friend mint operator+(const mint& lhs, const mint& rhs) {
        return mint(lhs) += rhs;
    }
    friend mint operator-(const mint& lhs, const mint& rhs) {
        return mint(lhs) -= rhs;
    }
    friend mint operator*(const mint& lhs, const mint& rhs) {
        return mint(lhs) *= rhs;
    }
    friend mint operator/(const mint& lhs, const mint& rhs) {
        return mint(lhs) /= rhs;
    }
    friend bool operator==(const mint& lhs, const mint& rhs) {
        return lhs._v == rhs._v;
    }
    friend bool operator!=(const mint& lhs, const mint& rhs) {
        return lhs._v != rhs._v;
    }

  private:
    unsigned int _v;
    static internal::barrett bt;
    static unsigned int umod() { return bt.umod(); }
};
template <int id> internal::barrett dynamic_modint<id>::bt(998244353);

using modint998244353 = static_modint<998244353>;
using modint1000000007 = static_modint<1000000007>;
using modint = dynamic_modint<-1>;

namespace internal {

template <class T>
using is_static_modint = std::is_base_of<internal::static_modint_base, T>;

template <class T>
using is_static_modint_t = std::enable_if_t<is_static_modint<T>::value>;

template <class> struct is_dynamic_modint : public std::false_type {};
template <int id>
struct is_dynamic_modint<dynamic_modint<id>> : public std::true_type {};

template <class T>
using is_dynamic_modint_t = std::enable_if_t<is_dynamic_modint<T>::value>;

}  // namespace internal

}  // namespace atcoder

#endif  // ATCODER_MODINT_HPP

#include <algorithm>

namespace atcoder {

template <int m, std::enable_if_t<(1 <= m)>* = nullptr>
struct combinatorics : static_modint<m> {
    using mint = static_modint<m>;

  public:

    combinatorics() = default;
    template <class T> combinatorics(T v) : mint(v) { }

    mint inv() const {
        assert(mint::val() != 0);
        if (mint::val() < preprocess_limit) {
            preprocess(mint::val());
            return _inv[mint::val()];
        }
        return mint::inv();
    }

    mint fact() const {
        preprocess(mint::val());
        return _fact[mint::val()];
    }

    mint fact_inv() const {
        preprocess(mint::val());
        return _ifact[mint::val()];
    }

    // Explicitly needed to create those static functions 
    // because the struct does not allow negative data 
    // processing (the value is stored in _v (unsigned int))

    static mint inv(int x) {
        if (x < 0) return 0;
        return combinatorics(x).inv();        
    }

    static mint fact(int x) {
        if (x < 0) return 0;
        return combinatorics(x).fact();
    }

    static mint fact_inv(int x) {
        if (x < 0) return 0;
        return combinatorics(x).fact_inv();
    }

    static mint P(int n, int r) {
        return fact(n) * fact_inv(n - r);
    }

    static mint C(int n, int r) {
        return P(n, r) * fact_inv(r);
    }
    
    // sqrt algorithm by Tonelli and Shanks
    // Reference: https://en.wikipedia.org/wiki/Tonelli%E2%80%93Shanks_algorithm
    // Time complexity: O(log^2 m)
    mint sqrt() const {
        // assert(prime);
        if (mint::val() == 0) return 0;
        if (m == 2) return mint::val();
        int e = 0, md = m - 1;
        if (mint::pow(md / 2) != 1) return 0;
        mint b = 1;
        while (b.mint::pow(md / 2) == 1) b += 1;
        while (md % 2 == 0) md >>= 1, e++;
        mint x = mint::pow((md - 1) / 2), y = (*this) * x * x;
        x *= (*this);
        mint z = b.mint::pow(md);
        while (y != 1) {
            int j = 0;
            mint t = y;
            while (t != 1) t *= t, j++;
            z = z.mint::pow(1LL << (e - j - 1));
            x *= z, z *= z, y *= z;
            e = j;
        }
        return mint(std::min(x.val(), m - x.val()));
    }

  private:
    static std::vector<unsigned int> _fact, _ifact, _inv;
    static unsigned int processed_until;
    static constexpr unsigned int preprocess_limit = (1 << 20);

    static void preprocess(int val) {
        if (val <= processed_until) return;
        _ifact.resize(val + 1);
        _fact.resize(val + 1);
        _inv.resize(val + 1);
        if (val == 0) return;
        for (long long i = processed_until + 1; i <= val; i++) {
            _fact[i] = _fact[i - 1] * i % m;
        }
        _ifact[val] = mint(_fact[val]).inv().val();
        for (long long i = val - 1; i >= processed_until; i--) {
            _ifact[i] = _ifact[i + 1] * (i + 1) % m;
            _inv[i + 1] = _ifact[i + 1] * 1ll * _fact[i] % m;
        }
        processed_until = val;
    }
};

} // namespace atcoder

#endif  // ATCODER_COMBINATORICS_HPP

// DO NOT FORGET TO INITIALIZE THE STATIC VARIABLES

// using comb = atcoder::combinatorics<998244353>;

// template<> std::vector<unsigned int> comb::_fact = {1};
// template<> std::vector<unsigned int> comb::_ifact = {1};
// template<> std::vector<unsigned int> comb::_inv = {1};
// template<> unsigned int comb::processed_until = 0;


using comb = atcoder::combinatorics<998244353>;

template<> std::vector<unsigned int> comb::_fact = {1};
template<> std::vector<unsigned int> comb::_ifact = {1};
template<> std::vector<unsigned int> comb::_inv = {1};
template<> unsigned int comb::processed_until = 0;


signed main() {

    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    int n;
    std::cin >> n;
    atcoder::matrix<comb> A(n, n);

    for (int i = 0; i < n; i++) {
        for (int j = 0, x; j < n; j++) {
            std::cin >> x;
            A[i][j] = x;
        }
    }

    auto inv = A.inverse();
    if (inv.rows_count == 0) {
        std::cout << -1 << std::endl;
    }
    else {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                std::cout << inv[i][j].val() << ' ';
            }
            std::cout << std::endl;
        }
    }
}