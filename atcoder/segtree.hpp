#ifndef ATCODER_SEGTREE_HPP
#define ATCODER_SEGTREE_HPP 1

#include <algorithm>
#include <cassert>
#include <functional>
#include <vector>
#include <bit>

namespace atcoder {

#ifdef COMMENT

int op(int a, int b) { return a + b; }
int e() { return 0; }

#endif

template <class S, S (*op)(S, S), S (*e)()> struct segtree {

  public:
    segtree() : segtree(0) {}
    explicit segtree(unsigned int N) : segtree(std::vector<S>(N, e())) {}
    explicit segtree(const std::vector<S> &v) : size(v.size()), n(std::bit_ceil(v.size())), d(2 * n, e()) {
        std::copy(v.begin(), v.end(), d.begin() + n);
        for (unsigned int i = n - 1; i; i--) {
            pull(i);
        }
    }

    void point_update(unsigned int idx, const S &x) {
        assert(0 <= idx && idx < size);
        for (d[idx += n] = x; idx /= 2;) {
            pull(idx);
        }
    }

    S range_query(unsigned int l, unsigned int r) {
        assert(0 <= l && l <= r && r < size);
        S res = e();
        for (l += n, r += n + 1; l < r; l /= 2, r /= 2) {
            if (l & 1) res = op(res, d[l++]);
            if (r & 1) res = op(res, d[--r]);
        }
        return res;
    }

    // returns maximum r >= l such that f(op(a[l], a[l+1], ..., a[r])) is true
    unsigned int max_right(unsigned int l, std::function<bool(S)> f) {
        assert(0 <= l && l < size);
        l += n;
        if (!f(d[l])) return (unsigned int)(-1);
        S sm = e();
        do {
            while (l % 2 == 0) l /= 2;
            if (!f(op(sm, d[l]))) {
                while (l < n) {
                    l = 2 * l;
                    if (f(op(sm, d[l]))) {
                        sm = op(sm, d[l++]);
                    }
                }
                return l - n - 1;
            }
            sm = op(sm, d[l++]);
        } while (l & (l - 1));
        return size - 1;
    }

    // returns minimum l <= r such that f(op(a[l], a[l+1], ..., a[r])) is true
    unsigned int min_left(unsigned int r, std::function<bool(S)> f) {
        assert(0 <= r && r < size);
        r += n;
        if (!f(d[r])) return (unsigned int)(-1);
        S sm = e();
        do {
            while (r > 1 && (r % 2)) r /= 2;
            if (!f(op(d[r], sm))) {
                while (r < n) {
                    r = 2 * r + 1;
                    if (f(op(d[r], sm))) {
                        sm = op(d[r--], sm);
                    }
                }
                return r - n + 1;
            }
            sm = op(d[r--], sm);
        } while (r & (r + 1));
        return 0;
    }

  private:

    unsigned int size, n;
    std::vector<S> d;

    void pull(unsigned int node) {
        d[node] = op(d[2 * node], d[2 * node + 1]);
    }
};

}  // namespace atcoder

#endif  // ATCODER_SEGTREE_HPP
