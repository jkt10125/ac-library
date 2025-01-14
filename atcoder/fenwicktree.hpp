#ifndef ATCODER_FENWICKTREE_HPP
#define ATCODER_FENWICKTREE_HPP 1

#include <cassert>
#include <vector>
#include <array>

namespace atcoder {

// Reference: https://en.wikipedia.org/wiki/Fenwick_tree
struct fenwick_tree {
    using T = long long;
  public:
    fenwick_tree() = default;
    explicit fenwick_tree(int n) : d(n + 1, {T(), T()}) {}

    void range_add(int l, int r, T x) {
        assert(0 <= l && l <= r && r < int(d.size()));
        set(l, x), set(r + 1, -x);
    }

    T range_sum(int l, int r) const {
        assert(0 <= l && l <= r && r < int(d.size()));
        return (get(r + 1) - get(l));
    }

  private:
    std::vector<std::array<T, 2>> d;

    void set(int i, T x) {
        for (int j = i + 1; j <= int(d.size()); j += j & -j) {
            d[j - 1][0] += x;
            d[j - 1][1] += x * i;
        }
    }

    T get(int i) const {
        T s0{}, s1{};
        for (int j = i + 1; j; j -= j & -j) {
            s0 += d[j - 1][0];
            s1 += d[j - 1][1];
        }
        return s0 * i - s1;
    }
};

}  // namespace atcoder

#endif  // ATCODER_FENWICKTREE_HPP
