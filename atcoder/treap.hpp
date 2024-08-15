#ifndef ATCODER_TREAP_HPP
#define ATCODER_TREAP_HPP

#include <vector>
#include <iostream>

namespace atcoder {

template <class T>
struct treap {

  public:
    treap() : root(new_node({}, 0)) {
        d[root].cnt = 0;
    }

    void insert(unsigned int idx, T val, unsigned int priority) {
        unsigned int l, r;
        split(root, l, r, idx);
        root = merge(merge(l, new_node(val, priority)), r);
    }

    void append(T val, unsigned int priority) {
        root = merge(root, new_node(val, priority));
    }

    void erase(unsigned int l, unsigned int r) {
        unsigned int a, b, c;
        split(root, b, c, r + 1);
        split(b, a, b, l);
        mark_for_deletion(b);
        root = merge(a, c);
    }

    void reverse(unsigned int l, unsigned int r) {
        unsigned int a, b, c;
        split(root, b, c, r + 1);
        split(b, a, b, l);
        d[b].rev ^= true;
        root = merge(merge(a, b), c);
    }

    void cyclic_shift(unsigned int l, unsigned int r, unsigned int k) {
        unsigned int a, b1, b2, c;
        split(root, b2, c, r + 1);
        split(b2, a, b1, l);
        unsigned int len = d[c].cnt;
        k = (k % len + len) % len;
        split(b1, b1, b2, len - k);
        root = merge(merge(a, merge(b2, b1)), c);
    }

    T& operator [] (unsigned int idx) {
        unsigned int a, b, c;
        split(root, b, c, idx + 1);
        split(b, a, b, idx);
        T &ret = d[b].val;
        root = merge(merge(a, b), c);
        return ret;
    }

    unsigned int size() {
        return d[root].cnt;
    }

    friend std::ostream& operator << (std::ostream& os, treap& t) {
        t.print_rec(os, t.root);
        return os;
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
        if (!d[id].rev) return;
        std::swap(d[id].left, d[id].right);
        if (d[id].left) d[d[id].left].rev ^= true;
        if (d[id].right) d[d[id].right].rev ^= true;
        d[id].rev = false;
    }

    void update_cnt(unsigned int id) {
        if (!id) return;
        d[id].cnt = 1 + d[d[id].left].cnt + d[d[id].right].cnt;
    }

    void split(unsigned int id, unsigned int &l, unsigned int &r, unsigned int k) {
        if (!id) return void(l = r = 0);
        push(id);
        if (d[d[id].left].cnt < k) {
            split(d[id].right, d[id].right, r, k - d[d[id].left].cnt - 1);
            l = id;
        }
        else {
            split(d[id].left, l, d[id].left, k);
            r = id;
        }
        update_cnt(id);
    }

    int merge(unsigned int l, unsigned int r) {
        if (!l || !r) return l + r;
        push(l);
        push(r);
        if (d[l].pri > d[r].pri) {
            d[l].right = merge(d[l].right, r);
            return update_cnt(l), l;
        } else {
            d[r].left = merge(l, d[r].left);
            return update_cnt(r), r;
        }
    }

    

    void print_rec(std::ostream& os, unsigned int id) {
        if (!id) return;
        print_rec(os, d[id].left);
        os << d[id].val << ' ';
        print_rec(os, d[id].right);
    }
};


}   // namespace atcoder

#endif  // ATCODER_TREAP_HPP