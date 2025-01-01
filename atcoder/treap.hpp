#ifndef ATCODER_TREAP_HPP
#define ATCODER_TREAP_HPP

#include <vector>
#include <cassert>
#include <iostream>

namespace atcoder {

template <typename T>
struct treap {
  public:

    treap() : d(1, {T(), 0, 0, 0, 0, false}), root(0) { }

    void insert(unsigned int idx, T val, unsigned int priority) {
        assert(0 <= idx && idx <= size());
        unsigned int l, r;
        split(root, l, r, idx);
        root = merge(merge(l, new_node(val, priority)), r);
    }

    void append(T val, unsigned int priority) {
        root = merge(root, new_node(val, priority));
    }

    void erase(unsigned int l, unsigned int r) {
        assert(0 <= l && l <= r && r < size());
        unsigned int a, b, c;
        split(root, b, c, r + 1);
        split(b, a, b, l);
        mark_for_deletion(b);
        root = merge(a, c);
    }

    void reverse(unsigned int l, unsigned int r) {
        assert(0 <= l && l <= r && r < size());
        unsigned int a, b, c;
        split(root, b, c, r + 1);
        split(b, a, b, l);
        d[b].rev ^= true;
        root = merge(merge(a, b), c);
    }

    void cyclic_shift(unsigned int l, unsigned int r, unsigned int k) {
        assert(0 <= l && l <= r && r < size());
        k %= (r - l + 1);
        unsigned int a, b1, b2, c;
        split(root, b2, c, r + 1);
        split(b2, a, root, l);
        split(root, b1, b2, r - l + 1 - k);
        root = merge(merge(a, merge(b2, b1)), c);
    }

    T &operator [] (unsigned int idx) {
        assert(0 <= idx && idx < size());
        unsigned int l, m, r;
        split(root, m, r, idx + 1);
        split(m, l, root, idx);
        m = root;
        root = merge(merge(l, root), r);
        return d[m].val;
    }

    unsigned int size() {
        return d[root].cnt;
    }

  private:
    
    struct node {
        T val;
        unsigned int cnt, pri, left, right;
        bool rev;
    };

    std::vector<node> d;
    std::vector<unsigned int> del;
    unsigned int root;

    unsigned int new_node(T val, unsigned int priority) {
        unsigned int id;
        if (del.empty()) {
            id = d.size();
            d.push_back({val, 1, priority, 0, 0, false});
        } else {
            id = del.back();
            del.pop_back();
            d[id] = {val, 1, priority, 0, 0, false};
        }
        return id;
    }

    void mark_for_deletion(unsigned int id) {
        if (!id) return;
        mark_for_deletion(d[id].left);
        mark_for_deletion(d[id].right);
        del.push_back(id);
    }

    void push(unsigned int id) {
        if (d[id].rev) {
            d[id].rev = false;
            std::swap(d[id].left, d[id].right);
            if (d[id].left) d[d[id].left].rev ^= 1;
            if (d[id].right) d[d[id].right].rev ^= 1;
        }
    }

    void update_node(unsigned int id) {
        if (!id) return;
        d[id].cnt = 1 + d[d[id].left].cnt + d[d[id].right].cnt;
    }

    void split(unsigned int id, unsigned int &l, unsigned int &r, unsigned int k) {
        if (!id) return void(l = r = 0);
        push(id);
        if (d[d[id].left].cnt < k) {
            split(d[id].right, d[id].right, r, k - d[d[id].left].cnt - 1);
            l = id;
        } else {
            split(d[id].left, l, d[id].left, k);
            r = id;
        }
        update_node(id);
    }

    int merge(unsigned int l, unsigned int r) {
        if (!l || !r) return l + r;
        push(l); push(r);
        if (d[l].pri > d[r].pri) {
            d[l].right = merge(d[l].right, r);
            return update_node(l), l;
        } else {
            d[r].left = merge(l, d[r].left);
            return update_node(r), r;
        }
    }
};

}   // namespace atcoder

#endif  // ATCODER_TREAP_HPP
