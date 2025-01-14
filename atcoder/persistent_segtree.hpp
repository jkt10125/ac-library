#ifndef ATCODER_PERSISTENT_SEGTREE_HPP
#define ATCODER_PERSISTENT_SEGTREE_HPP 1

// #include "internal_bit.hpp"

#include <vector>
#include <cassert>
#include <iostream>
#include <bit>

namespace atcoder {

#ifdef COMMENT

int op(int a, int b) { return a + b; }
int e() { return 0; }

#endif

template <class S, S (*op)(S, S), S (*e)()> 
struct persistent_segtree {
  public:

    persistent_segtree() : persistent_segtree(0) {}
    explicit persistent_segtree(unsigned int N) : persistent_segtree(std::vector<S>(N, e())) {}
    explicit persistent_segtree(const std::vector<S> &v) : size(v.size()), 
                                                           n(std::bit_ceil(v.size())), 
                                                           d(2 * n, e()), 
                                                           lch(2 * n), 
                                                           rch(2 * n), 
                                                           root(1, 1) {
        std::copy(v.begin(), v.end(), d.begin() + n);
        for (unsigned int i = n - 1; i; i--) {
            lch[i] = 2 * i;
            rch[i] = 2 * i + 1;
            pull(i);
        }
    }

    void point_update(unsigned int ver, unsigned int idx, const S &x) {
        assert(0 <= idx && idx < size);
        assert(0 <= ver && ver < root.size());
        root.push_back(update(idx, x, root[ver], 0, n - 1));
    }

    S range_query(unsigned int ver, unsigned int l, unsigned int r) {
        assert(0 <= l && l <= r && r < size);
        assert(0 <= ver && ver < root.size());
        return query(l, r, root[ver], 0, n - 1);
    }

    void print(unsigned int ver, std::ostream &os) {
        for (unsigned int i = 0; i < size; i++) {
            os << range_query(ver, i, i) << ' ';
        }
    }

  private:
    unsigned int size, n;
    std::vector<S> d;
    std::vector<unsigned int> lch, rch;
    std::vector<unsigned int> root;

    unsigned int new_node() {
        d.push_back(e());
        lch.push_back(0);
        rch.push_back(0);
        return (d.size() - 1);
    }

    unsigned int update(unsigned int idx, const S &x,
                        unsigned int node, unsigned int l, unsigned int r) {
        unsigned int id = new_node();
        if (l == r) {
            d[id] = x;
        } else {
            unsigned int m = (l + r) / 2;
            if (m < idx) {
                lch[id] = lch[node];
                rch[id] = update(idx, x, rch[node], m + 1, r);
            } else {
                rch[id] = rch[node];
                lch[id] = update(idx, x, lch[node], l, m);
            }
            pull(id);
        }
        return id;
    }

    S query(unsigned int L, unsigned int R, 
            unsigned int node, unsigned int l, unsigned int r) {
        if (r < L || R < l) {
            return e();
        } else if (L <= l && r <= R) {
            return d[node];
        } else {
            unsigned int m = (l + r) / 2;
            return op(query(L, R, lch[node], l, m), query(L, R, rch[node], m + 1, r));
        }
    }

    void pull(unsigned int node) {
        d[node] = op(d[lch[node]], d[rch[node]]);
    }
};

}

#endif