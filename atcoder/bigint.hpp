#ifndef ATCODER_BIGINT_HPP
#define ATCODER_BIGINT_HPP 1

#include <vector>
#include <string>
#include <algorithm>

#include "atcoder/convolution"

namespace atcoder {

// base must be power of 10
struct bigint {
  public:
    bigint() : base(10), is_negative(false), num(1, 0) {
        calc_base_digits();
    }

    bigint(const bigint &x) = default;
    
    bigint& operator=(const bigint &x) = default;

    bigint(long long x) : base(10), is_negative(x < 0) {
        if (x < 0) x = -x;
        if (x == 0) num.assign(1, 0);
        else {
            num.clear();
            while (x) {
                num.push_back(x % base);
                x /= base;
            }
        }
        calc_base_digits();
    }
    
    bigint(const std::string &s) : base(10), is_negative(s[0] == '-') {
        int start = is_negative;
        num.clear();
        for (int i = int(s.size()) - 1; i >= start; i--) {
            num.push_back(s[i] - '0');
        }
        while (num.size() > 1 && num.back() == 0) num.pop_back();
        calc_base_digits();
    }
    
    bigint(const std::vector<long long> &v, bool sign) : base(10), is_negative(sign), num(v) {
        while (num.size() > 1 && num.back() == 0) num.pop_back();
        calc_base_digits();
    }
    
    bool operator==(const bigint &x) const {
        return is_negative == x.is_negative && num == x.num;
    }
    
    bool operator<(const bigint &x) const {
        if (is_negative != x.is_negative) return is_negative;
        if (num.size() != x.num.size()) return (num.size() < x.num.size()) ^ is_negative;
        return std::lexicographical_compare(num.rbegin(), num.rend(), x.num.rbegin(), x.num.rend()) ^ is_negative;
    }
    
    bool abs_compare_lt(const bigint &x) const {
        if (num.size() != x.num.size()) return num.size() < x.num.size();
        return std::lexicographical_compare(num.rbegin(), num.rend(), x.num.rbegin(), x.num.rend());
    }
    
    void negate() { is_negative ^= 1; }
    
    bigint operator-() const {
        bigint ret = *this;
        ret.negate();
        return ret;
    }
    
    bigint& operator+=(const bigint &x) {
        // static_assert(base == x.base, "mismatched base");
        if (is_negative == x.is_negative) {
            if (num.size() < x.num.size()) num.resize(x.num.size(), 0);
            int carry = 0;
            for (size_t i = 0; i < x.num.size() || carry; i++) {
                if (i == num.size()) num.push_back(0);
                num[i] += carry + (i < x.num.size() ? x.num[i] : 0);
                carry = num[i] >= base;
                if (carry) num[i] -= base;
            }
        } else {
            if (!abs_compare_lt(x)) {
                int carry = 0;
                for (size_t i = 0; i < x.num.size() || carry; i++) {
                    num[i] -= carry + (i < x.num.size() ? x.num[i] : 0);
                    carry = num[i] < 0;
                    if (carry) num[i] += base;
                }
                while (num.size() > 1 && num.back() == 0) num.pop_back();
            } else {
                bigint tmp = x;
                tmp -= *this;
                *this = tmp;
                is_negative = x.is_negative;
            }
        }
        return *this;
    }
    
    bigint& operator-=(const bigint &x) {
        return *this += -x;
    }
    
    bigint& operator*=(int x) {
        // static_assert(0 <= x && x < base, "out of range");
        if (x == 0) return *this = 0;
        if (x < 0) is_negative ^= 1, x = -x;
        int carry = 0;
        for (size_t i = 0; i < num.size() || carry; i++) {
            if (i == num.size()) num.push_back(0);
            long long cur = num[i] * 1LL * x + carry;
            carry = cur / base;
            num[i] = cur % base;
        }
        while (num.size() > 1 && num.back() == 0) num.pop_back();
        return *this;
    }
    
    bigint& operator*=(const bigint &x) {
        return *this = *this * x;
    }
    
    bigint operator*(const bigint &x) const {
        // static_assert(base == x.base, "mismatched base");
        if (num.size() == 1 && num[0] == 0) return 0;
        if (x.num.size() == 1 && x.num[0] == 0) return 0;
        std::vector<long long> a(num.begin(), num.end());
        std::vector<long long> b(x.num.begin(), x.num.end());
        std::vector<long long> c = convolution_ll(std::move(a), std::move(b));
        int carry = 0;
        for (int i = 0; i < int(c.size()); i++) {
            c[i] += carry;
            if (c[i] < base) carry = 0;
            else carry = c[i] / base, c[i] %= base;
        }
        if (carry) c.push_back(carry);
        return bigint(c, is_negative ^ x.is_negative);
    }

    std::string to_string() const {
        std::string ans;
        if (is_negative) {
            ans.push_back('-');
        }
        for (int i = int(num.size()) - 1; i >= 0; --i) {
            std::string&& s = std::to_string(num[i]);
            if (i < int(num.size()) - 1) {
                s = std::string(base_digits - int(s.size()), '0') + s;
            }
            ans += s;
        }
        return ans;
    }

    // input output
    friend std::istream& operator>>(std::istream& is, bigint& x) {
        std::string s;
        is >> s;
        x = bigint(s);
        return is;
    }

    friend std::ostream& operator<<(std::ostream& os, bigint& x) {
        os << x.to_string();
        return os;
    }
    
  private:
    int base, base_digits;
    bool is_negative;
    std::vector<long long> num;

    void calc_base_digits() {
        base_digits = 0;
        for (int x = base - 1; x; x /= 10) base_digits++;
    }
};

} // namespace atcoder

#endif // ATCODER_BIGINT_HPP