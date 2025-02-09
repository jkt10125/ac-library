#ifndef ATCODER_LAZYSEGTREE_HPP
#define ATCODER_LAZYSEGTREE_HPP 1

#include <vector>
#include <functional>
#include <cassert>
#include <bit>

namespace atcoder {

#ifdef RANGE_ADD_RANGE_MIN

struct Tag {
    long long v;
    Tag() : v(0) { }
    Tag(long long x) : v(x) { }
    void apply(const Tag &t) {
        v += t.v;
    }
    bool operator == (const Tag &t) const {
        return v == t.v;
    }
};
 
struct Info {
    long long v;
    Info() : v(0) { }
    Info(long long x) : v(x) { }
    void apply(const Tag &t, unsigned int size) {
        v += t.v;
    }
    static Info merge(const Info &lhs, const Info &rhs) {
        return std::min(lhs.v, rhs.v);
    }
    static Info identity() {
        return 1000000000000000007ll;
    }
};

#endif

#ifdef RANGE_ADD_RANGE_SUM

struct Tag {
    long long v;
    Tag() : v(0) { }
    Tag(long long x) : v(x) { }
    void apply(const Tag &t) {
        v += t.v;
    }
    bool operator == (const Tag &t) const {
        return v == t.v;
    }
};
 
struct Info {
    long long v;
    Info() : v(0) { }
    Info(long long x) : v(x) { }
    void apply(const Tag &t, unsigned int size) {
        v += t.v * size;
    }
    static Info merge(const Info &lhs, const Info &rhs) {
        return Info(lhs.v + rhs.v);
    }
    static Info identity() {
        return 0ll;
    }
};

#endif

#ifdef RANGE_SET_RANGE_MIN

struct Tag {
    long long v;
    Tag() : v(1000000000000000007ll) { }
    Tag(long long x) : v(x) { }
    void apply(const Tag &t) {
        v = t.v;
    }
    bool operator == (const Tag &t) const {
        return v == t.v;
    }
};
 
struct Info {
    long long v;
    Info() : v(0) { }
    Info(long long x) : v(x) { }
    void apply(const Tag &t, unsigned int size) {
        v = t.v;
    }
    static Info merge(const Info &lhs, const Info &rhs) {
        return std::min(lhs.v, rhs.v);
    }
    static Info identity() {
        return 1000000000000000007ll;
    }
};

#endif

#ifdef RANGE_SET_RANGE_SUM

struct Tag {
    long long v;
    Tag() : v(1000000000000000007ll) { }
    Tag(long long x) : v(x) { }
    void apply(const Tag &t) {
        v = t.v;
    }
    bool operator == (const Tag &rhs) const {
        return v == rhs.v;
    }
};
 
struct Info {
    long long v;
    Info() : v(0) { }
    Info(long long x) : v(x) { }
    void apply(const Tag &t, unsigned int size) {
        v = t.v * size;
    }
    static Info merge(const Info &lhs, const Info &rhs) {
        return lhs.v + rhs.v;
    }
    static Info identity() {
        return 0ll;
    }
};

#endif

#ifdef COMMENT

struct Tag {
    Tag() { }
    Tag(int x) { }
    void apply(const Tag &t) { }
    bool operator == (const Tag &t) { }
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

    lazysegtree(unsigned int N) : n(std::bit_ceil(N)), tree(2 * n), lazy(2 * n) { }

    void point_update(unsigned int idx, const Info &x) {
        assert(0 <= idx && idx < n);
        point_update(idx, x, 1, 0, n - 1);
    }
    
    void range_update(unsigned int L, unsigned int R, const Tag &t) {
        assert(0 <= L && L <= R && R < n);
        range_update(L, R, t, 1, 0, n - 1);
    }

    Info range_query(unsigned int L, unsigned int R) {
        assert(0 <= L && L <= R && R < n);
        return range_query(L, R, 1, 0, n - 1);
    }

    unsigned int min_left(unsigned int L, unsigned int R, std::function<bool(Info)> pred) {
        assert(0 <= L && L <= R && R < n);
        return min_left(L, R, pred, 1, 0, n - 1);
    }

    unsigned int max_right(unsigned int L, unsigned int R, std::function<bool(Info)> pred) {
        assert(0 <= L && L <= R && R < n);
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