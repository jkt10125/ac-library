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

    matrix(std::vector<std::vector<T>> d_) : d(std::move(d_)), rows_count((*this).d.size()), columns_count((*this).d[0].size()) { }

    std::vector<T>& operator[](unsigned int i) { return d[i]; }
    const std::vector<T>& operator[](unsigned int i) const { return d[i]; }

    static matrix identity(unsigned int n) {
        matrix mat(n, n);
        for (int i = 0; i < int(n); i++) mat[i][i] = 1;
        return mat;
    }

    bool empty() const { return rows_count == 0 || columns_count == 0; }

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

  private:
    std::vector<std::vector<T>> d;
    unsigned int rows_count, columns_count;
};

}  // namespace atcoder

#endif  // ATCODER_MATRIX_HPP