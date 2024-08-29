#ifndef ATCODER_LAZYSEGTREE_HPP
#define ATCODER_LAZYSEGTREE_HPP 1

#include "internal_bit.hpp"

#include <vector>
#include <type_traits>

namespace atcoder {

#ifdef COMMENT

struct Tag {
    Tag() { }
    Tag(int x) { }
    void apply(const Tag &t) { }
    bool operator == (const Tag &t) { }
    // Always return true in the == operator if you do not want lazy propagation
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

    void pointUpdate(int idx, const Info &x) {
        pointUpdate(idx, x, 1, 0, n - 1);
    }
    
    void rangeUpdate(int L, int R, const Tag &t) {
        rangeUpdate(L, R, t, 1, 0, n - 1);
    }

    Info rangeQuery(int L, int R) {
        return rangeQuery(L, R, 1, 0, n - 1);
    }

    template <typename F>
    int findFirst(int L, int R, F pred) {
        return findFirst(L, R, pred, 1, 0, n - 1);
    }

    template <typename F>
    int findLast(int L, int R, F pred) {
        return findLast(L, R, pred, 1, 0, n - 1);
    }

    int findKth(int L, int R, int k) {
        return findKth(L, R, k, 1, 0, n - 1);
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

    void pointUpdate(unsigned int idx, const Info &x, 
                     unsigned int node, unsigned int l, unsigned int r) {
        if (l == r) {
            tree[node] = x;
            return;
        }
        unsigned int m = (l + r) / 2;
        push(node, r - l + 1);
        if (m < idx) pointUpdate(idx, x, 2 * node + 1, m + 1, r);
        else pointUpdate(idx, x, 2 * node, l, m);
        pull(node);
    }

    void rangeUpdate(unsigned int L, unsigned int R, const Tag &t, 
                     unsigned int node, unsigned int l, unsigned int r) {
        if (r < L || R < l) return;
        if (L <= l && r <= R) {
            apply(node, t, r - l + 1);
            return;
        }
        unsigned int m = (l + r) / 2;
        push(node, r - l + 1);
        rangeUpdate(L, R, t, 2 * node, l, m);
        rangeUpdate(L, R, t, 2 * node + 1, m + 1, r);
        pull(node);
    }

    Info rangeQuery(unsigned int L, unsigned int R, 
                    unsigned int node, unsigned int l, unsigned int r) {
        if (r < L || R < l) return Info::identity();
        if (L <= l && r <= R) {
            return tree[node];
        }
        unsigned int m = (l + r) / 2;
        push(node, r - l + 1);
        return Info::merge(rangeQuery(L, R, 2 * node, l, m), rangeQuery(L, R, 2 * node + 1, m + 1, r));
    }

    template <typename F>
    unsigned int findFirst(unsigned int L, unsigned int R, F pred, 
                           unsigned int node, unsigned int l, unsigned int r) {
        if (r < L || R < l || !pred(node)) return unsigned(-1);
        if (l == r) return l;
        unsigned int m = (l + r) / 2;
        push(node, r - l + 1);
        unsigned int res = findFirst(L, R, pred, 2 * node, l, m);
        if (res == unsigned(-1)) {
            res = findFirst(L, R, pred, 2 * node + 1, m + 1, r);
        }
        return res;
    }

    template <typename F>
    unsigned int findLast(unsigned int L, unsigned int R, F pred, 
                          unsigned int node, unsigned int l, unsigned int r) {
        if (r < L || R < l || !pred(node)) return unsigned(-1);
        if (l == r) return l;
        unsigned int m = (l + r) / 2;
        push(node, r - l + 1);
        unsigned int res = findLast(L, R, pred, 2 * node + 1, m + 1, r);
        if (res == unsigned(-1)) {
            res = findLast(L, R, pred, 2 * node, l, m);
        }
        return res;
    }

    unsigned int findKth(unsigned int L, unsigned int R, unsigned int k, 
                         unsigned int node, unsigned int l, unsigned int r) {
        if (r < L || R < l || tree[node].cnt < k) return unsigned(-1);
        if (l == r) return l;
        int m = (l + r) / 2;
        push(node, r - l + 1);
        int res = findKth(L, R, k, 2 * node, l, m);
        if (res == unsigned(-1)) {
            res = findKth(L, R, k - tree[2 * node].cnt, 2 * node + 1, m + 1, r);
        }
        return res;
    }

    unsigned int n;
    std::vector<Info> tree;
    std::vector<Tag> lazy;
};

}   // namespace atcoder

#endif  // ATCODER_LAZYSEGTREE_HPP