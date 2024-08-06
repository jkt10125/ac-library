#ifndef ATCODER_COMBINATORICS_HPP
#define ATCODER_COMBINATORICS_HPP 1

#include <algorithm>

#include "modint.hpp"

namespace atcoder {

template <int m, std::enable_if_t<(1 <= m)>* = nullptr>
struct combinatorics : static_modint<m> {
    using mint = static_modint<m>;

  public:

    combinatorics() = default;
    template <class T> combinatorics(T v) : mint(v) { }

    mint inv() const {
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

#endif

// DO NOT FORGET TO INITIALIZE THE STATIC VARIABLES

using comb = atcoder::combinatorics<998244353>;

template<> std::vector<unsigned int> comb::_fact = {1};
template<> std::vector<unsigned int> comb::_ifact = {1};
template<> std::vector<unsigned int> comb::_inv = {1};
template<> unsigned int comb::processed_until = 0;
