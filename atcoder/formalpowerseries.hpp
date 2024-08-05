#ifndef ATCODER_FORMAL_POWER_SERIES_HPP
#define ATCODER_FORMAL_POWER_SERIES_HPP 1

// Reference: https://ei1333.github.io/luzhiled/snippets/math/formal-power-series.html

#include <vector>
#include <algorithm>

#include "atcoder/modint"
#include "atcoder/convolution"

namespace atcoder {

using T = modint998244353;
// template <class T>
struct formal_power_series : std::vector<T> {
    using std::vector<T>::vector;
    using P = formal_power_series;

    void fix() {
        while (this->size() && this->back() == T(0)) this->pop_back();
    }

    P operator+(const P& r) const { return P(*this) += r; }
    P operator-(const P& r) const { return P(*this) -= r; }
    P operator*(const P& r) const { return P(*this) *= r; }
    P operator/(const P& r) const { return P(*this) /= r; }
    P operator%(const P& r) const { return P(*this) %= r; }

    P operator+(const T& r) const { return P(*this) += r; }
    P operator-(const T& r) const { return P(*this) -= r; }
    P operator*(const T& r) const { return P(*this) *= r; }
    P operator/(const T& r) const { return P(*this) /= r; }

    P& operator+=(const P& r) {
        if (this->size() < r.size()) this->resize(r.size());
        for (int i = 0; i < int(r.size()); ++i) (*this)[i] += r[i];
        return *this;
    }

    P& operator-=(const P& r) {
        if (this->size() < r.size()) this->resize(r.size());
        for (int i = 0; i < int(r.size()); ++i) (*this)[i] -= r[i];
        return *this;
    }

    P& operator*=(const P& r) {
        auto&& tmp = convolution(*this, r);
        *this = P(tmp.begin(), tmp.end());
        return *this;
    }

    P& operator%=(const P& r) {
        *this -= *this / r * r;
        return *this;
    }

    P& operator/=(const P& r) {
        if (this->size() < r.size()) {
            this->clear();
        } else {
            int n = int(this->size() - r.size()) + 1;
            *this = (rev().pre(n) * r.rev().inv(n)).pre(n).rev(n);
        }
        return *this;
    }


    P& operator+=(const T& r) {
        if (this->empty()) this->resize(1);
        (*this)[0] += r;
        return *this;
    }

    P& operator-=(const T& r) {
        if (this->empty()) this->resize(1);
        (*this)[0] -= r;
        return *this;
    }

    P& operator*=(const T& r) {
        for (T &x : *this) x *= r;
        return *this;
    }

    P& operator/=(const T& r) {
        for (T &x : *this) x /= r;
        return *this;
    }

    P operator-() const {
        P ret(*this);
        for (T &x : ret) x = -x;
        return ret;
    }

    P operator>>(int sz) const {
        if (int(this->size()) <= sz) return P();
        return P(this->begin() + sz, this->end());
    }

    P operator<<(int sz) const {
        if (this->empty()) return P();
        P ret(*this);
        ret.insert(ret.begin(), sz, T(0));
        return ret;
    }

    P differential() const {
        int n = int(this->size());
        P ret(std::max(0, n - 1));
        for (int i = 1; i < n; ++i) {
            ret[i - 1] = (*this)[i] * T(i);
        }
        return ret;
    }

    P integral() const {
        int n = int(this->size());
        P ret(n + 1, T(0));
        for (int i = 0; i < n; ++i) {
            ret[i + 1] = (*this)[i] / T(i + 1);
        }
        return ret;
    }

    P rev(int sz = -1) const {
        P ret(*this);
        if (sz == -1) ret.resize(sz, T(0));
        std::reverse(ret.begin(), ret.end());
        return ret;
    }

    P pre(int sz) const {
        return P(this->begin(), this->begin() + std::min(int(this->size()), sz));
    }

    P inv(int deg) const {
        if (deg == 0) return P();
        assert(!this->empty() && (*this)[0] != T(0));
        if (deg == -1) deg = int(this->size());
        P ret({T(1) / (*this)[0]});
        for (int i = 1; i < deg; i *= 2) {
            auto&& h = (pre(2 * i) * ret).pre(2 * i) >> i;
            auto&& tmp = (-h * ret).pre(i);
            ret.insert(ret.end(), tmp.begin(), tmp.end());
            ret.resize(2 * i);
        }
        return ret.pre(deg);
    }

    P log(int len = -1) const {
        if (len == 0) return P();
        assert(!this->empty() && (*this)[0] == T(1));
        if (len == -1) len = int(this->size());
        return (differential() * inv(len)).pre(len - 1).integral();
    }

    P sqrt(int deg = -1) const {
        if (this->empty()) return P();
        if (deg == -1) deg = int(this->size());
        if ((*this)[0] == T(0)) {
            for (int i = 1; i < int(this->size()); ++i) {
                if ((*this)[i] != T(0)) {
                    if (i % 2 == 1 || deg - i / 2 <= 0) return P();
                    return (*this >> i).sqrt(deg - i / 2) << (i / 2);
                }
            }
            return P();
        }
        T sqrtf0 = (*this)[0].sqrt();
        if (sqrtf0 == T(0)) return P();

        P y = (*this) / (*this)[0], ret({T(1)});
        T inv2 = T(1) / T(2);
        for (int i = 1; i < deg; i *= 2) {
            ret = (ret + y.pre(2 * i) * ret.inv(2 * i)) * inv2;
        }
        return ret.pre(deg) * sqrtf0;
    }

    P exp(int deg = -1) const {
        assert(this->empty() || (*this)[0] == T(0));
        if (deg == -1) deg = int(this->size());
        P ret({T(1)});
        for (int i = 1; i < deg; i *= 2) {
            ret = (ret * (pre(2 * i) + T(1) - ret.log(2 * i))).pre(2 * i);
        }
        return ret.pre(deg);
    }

    P pow(long long k, int deg = -1) const {
        const int n = int(this->size());
        if (deg == -1) deg = n;
        if (k == 0) {
            P ret(deg, T(0));
            if (deg >= 1) ret[0] = T(1);
            return ret;
        }
        for (int i = 0; i < n; ++i) {
            if ((*this)[i] != T(0)) {
                T rev = T(1) / (*this)[i];
                P C = (*this) * rev, D(n - i);
                for (int j = i; j < n; j++) {
                    D[j - i] = C[j];
                }
                D = (D.log(deg) * T(k)).exp(deg) * (*this)[i].pow(k);
                if (__int128(k) * i > deg) return P();
                P E(deg);
                long long S = i * k;
                for (int j = 0; j + S < deg && j < (int)D.size(); ++j) {
                    E[j + S] = D[j];
                }
                E.fix();
                return E;
            }
        }
        return *this;
    }

    P shift(T c) const {
        const int n = int(this->size());
        P ret = * this;
        for (int i = 0; i < n; i++) {
            ret[i] *= T(i).fact();
        }
        std::reverse(ret.begin(), ret.end());
        P exp_cx(n, 1);
        for (int i = 1; i < n; ++i) {
            exp_cx[i] = exp_cx[i - 1] * c * T(i).inv();
        }
        ret = ret * exp_cx;
        ret.resize(n);
        std::reverse(ret.begin(), ret.end());
        for (int i = 0; i < n; ++i) {
            ret[i] *= T(i).fact_inv();
        }
        return ret;
    }

    T eval(T x) const {
        T ret = 0, w = 1;
        for (auto& v : *this) {
            ret += w * v;
            w *= x;
        }
        return ret;
    }
};

} // namespace atcoder

#endif