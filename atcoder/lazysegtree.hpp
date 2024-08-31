#ifndef ATCODER_LAZYSEGTREE_HPP
#define ATCODER_LAZYSEGTREE_HPP 1

#include "internal_bit.hpp"

#include <vector>
#include <functional>
#include <type_traits>

namespace atcoder {

#ifdef COMMENT

struct Tag {
    Tag() { }
    Tag(int x) { }
    void apply(const Tag &t) { }
    bool operator == (const Tag &t) { }
    // return true in the == operator if you do not want lazy propagation
};

struct Info {
    Info() { }
    Info(int x) { }
    void apply(const Tag &t, unsigned int size) { }
    static Info merge(const Info &lhs, const Info &rhs) { }
    static Info identity() { }
};

#endif

template <typename Info, typename Tag>
struct lazysegtree {
  public:

    lazysegtree(unsigned int N) : n(internal::bit_ceil(N)), tree(2 * n), lazy(2 * n) { }

    void point_update(int idx, const Info &x) {
        point_update(idx, x, 1, 0, n - 1);
    }
    
    void range_update(int L, int R, const Tag &t) {
        range_update(L, R, t, 1, 0, n - 1);
    }

    Info range_query(int L, int R) {
        return range_query(L, R, 1, 0, n - 1);
    }

    unsigned int min_left(int L, int R, std::function<bool(Info)> pred) {
        return min_left(L, R, pred, 1, 0, n - 1);
    }

    unsigned int max_right(int L, int R, std::function<bool(Info)> pred) {
        return max_right(L, R, pred, 1, 0, n - 1);
    }

  private:
    
    void apply(unsigned int node, const Tag &t, unsigned int size) {
        tree[node].apply(t, size);
        lazy[node].apply(t);
    }

    void push(unsigned int node, unsigned int size) {
        if (lazy[node] == Tag()) return;
        apply(2 * node, lazy[node], size / 2);
        apply(2 * node + 1, lazy[node], size / 2);
        lazy[node] = Tag();
    }

    void pull(unsigned int node) {
        tree[node] = Info::merge(tree[2 * node], tree[2 * node + 1]);
    }

    void point_update(unsigned int idx, const Info &x, 
                     unsigned int node, unsigned int l, unsigned int r) {
        if (l == r) {
            tree[node] = x;
            return;
        }
        unsigned int m = (l + r) / 2;
        push(node, r - l + 1);
        if (m < idx) point_update(idx, x, 2 * node + 1, m + 1, r);
        else point_update(idx, x, 2 * node, l, m);
        pull(node);
    }

    void range_update(unsigned int L, unsigned int R, const Tag &t, 
                     unsigned int node, unsigned int l, unsigned int r) {
        if (r < L || R < l) return;
        if (L <= l && r <= R) {
            apply(node, t, r - l + 1);
            return;
        }
        unsigned int m = (l + r) / 2;
        push(node, r - l + 1);
        range_update(L, R, t, 2 * node, l, m);
        range_update(L, R, t, 2 * node + 1, m + 1, r);
        pull(node);
    }

    Info range_query(unsigned int L, unsigned int R, 
                    unsigned int node, unsigned int l, unsigned int r) {
        if (r < L || R < l) return Info::identity();
        if (L <= l && r <= R) {
            return tree[node];
        }
        unsigned int m = (l + r) / 2;
        push(node, r - l + 1);
        return Info::merge(range_query(L, R, 2 * node, l, m), range_query(L, R, 2 * node + 1, m + 1, r));
    }

    unsigned int min_left(unsigned int L, unsigned int R, std::function<bool(Info)> pred, 
                           unsigned int node, unsigned int l, unsigned int r) {
        if (r < L || R < l || !pred(tree[node])) return (unsigned int)(-1);
        if (l == r) return l;
        unsigned int m = (l + r) / 2;
        push(node, r - l + 1);
        unsigned int res = min_left(L, R, pred, 2 * node, l, m);
        if (res == (unsigned int)(-1)) {
            res = min_left(L, R, pred, 2 * node + 1, m + 1, r);
        }
        return res;
    }

    unsigned int max_right(unsigned int L, unsigned int R, std::function<bool(Info)> pred, 
                          unsigned int node, unsigned int l, unsigned int r) {
        if (r < L || R < l || !pred(tree[node])) return (unsigned int)(-1);
        if (l == r) return l;
        unsigned int m = (l + r) / 2;
        push(node, r - l + 1);
        unsigned int res = max_right(L, R, pred, 2 * node + 1, m + 1, r);
        if (res == (unsigned int)(-1)) {
            res = max_right(L, R, pred, 2 * node, l, m);
        }
        return res;
    }

    unsigned int n;
    std::vector<Info> tree;
    std::vector<Tag> lazy;
};

}   // namespace atcoder

#endif  // ATCODER_LAZYSEGTREE_HPP