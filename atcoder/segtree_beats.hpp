#ifndef ATCODER_SEGTREE_BEATS_HPP
#define ATCODER_SEGTREE_BEATS_HPP 1

// #include "internal_bit.hpp"

#include <vector>
#include <array>
#include <climits>
#include <cassert>
#include <bit>

namespace atcoder {

struct segtree_beats {
  public:
    
    segtree_beats() : segtree_beats(0) {}
    explicit segtree_beats(unsigned int N) : segtree_beats(std::vector<int>(N, INT_MIN)) {}
    explicit segtree_beats(const std::vector<int>& v) : size(v.size()), 
                                                        n(std::bit_ceil(v.size())), 
                                                        d(2 * n), 
                                                        lz(2 * n, INT_MAX) {
        for (unsigned int i = n; i < 2 * n; i++) {
            d[i] = Info(v[i - n]);
        }
        for (unsigned int i = n - 1; i; i--) {
            pull(i);
        }
    }

    void point_update(unsigned int idx, int x) {
        assert(0 <= idx && idx < size);
        point_update(idx, x, 1, 0, n - 1);
    }

    void range_update(unsigned int L, unsigned int R, int x) {
        assert(0 <= L && L <= R && R < size);
        range_update(L, R, x, 1, 0, n - 1);
    }

    long long range_sum(unsigned int L, unsigned int R) {
        assert(0 <= L && L <= R && R < size);
        return range_query(L, R, 1, 0, n - 1).sum;
    }

  private:

    struct Info {
        int max1, max2;
        long long sum;
        unsigned int max1_cnt;

        Info() : max1(INT_MIN), max2(INT_MIN), sum(0), max1_cnt(0) {}
        Info(long long v) : max1(v), max2(INT_MIN), sum(v), max1_cnt(1) {}

        static Info merge(const Info& a, const Info& b) {
            Info res;
            res.sum = a.sum + b.sum;
            if (a.max1 > b.max1) {
                res.max1 = a.max1;
                res.max1_cnt = a.max1_cnt;
                res.max2 = std::max(a.max2, b.max1);
            } else if (a.max1 < b.max1) {
                res.max1 = b.max1;
                res.max1_cnt = b.max1_cnt;
                res.max2 = std::max(a.max1, b.max2);
            } else {
                res.max1 = a.max1;
                res.max1_cnt = a.max1_cnt + b.max1_cnt;
                res.max2 = std::max(a.max2, b.max2);
            }
            return res;
        }

        void apply(int x) {
            if (max2 < x && x < max1) {
                sum += (long long)(x - max1) * max1_cnt;
                max1 = x;
            }
        }
    };

    void apply(unsigned int node, int x) {
        d[node].apply(x);
        lz[node] = std::min(lz[node], x);
    }

    void push(unsigned int node) {
        if (lz[node] == INT_MAX) return;
        apply(2 * node, lz[node]);
        apply(2 * node + 1, lz[node]);
        lz[node] = INT_MAX;
    }

    void pull(unsigned int node) {
        d[node] = Info::merge(d[2 * node], d[2 * node + 1]);
    }

    void point_update(unsigned int idx, int x, 
                      unsigned int node, unsigned int l, unsigned int r) {
        if (l == r) {
            d[node] = Info(x);
            return;
        }
        push(node);
        unsigned int m = (l + r) / 2;
        if (idx <= m) point_update(idx, x, 2 * node, l, m);
        else point_update(idx, x, 2 * node + 1, m + 1, r);
        pull(node);
    }

    void range_update(unsigned int L, unsigned int R, int x, 
                      unsigned int node, unsigned int l, unsigned int r) {
        if (R < l || r < L || d[node].max1 <= x) return;
        if (L <= l && r <= R && d[node].max2 < x) {
            apply(node, x);
            return;
        }
        push(node);
        unsigned int m = (l + r) / 2;
        range_update(L, R, x, 2 * node, l, m);
        range_update(L, R, x, 2 * node + 1, m + 1, r);
        pull(node);
    }

    Info range_query(unsigned int L, unsigned int R, 
                     unsigned int node, unsigned int l, unsigned int r) {
        if (R < l || r < L) return Info();
        if (L <= l && r <= R) return d[node];
        push(node);
        unsigned int m = (l + r) / 2;
        return Info::merge(range_query(L, R, 2 * node, l, m), range_query(L, R, 2 * node + 1, m + 1, r));
    }

    unsigned int size, n;
    std::vector<Info> d;
    std::vector<int> lz;
};

}   // namespace atcoder

#endif  // ATCODER_SEGTREE_BEATS_HPP