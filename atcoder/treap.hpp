#ifndef ATCODER_TREAP_HPP
#define ATCODER_TREAP_HPP

#include <vector>
#include <iostream>

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

namespace atcoder {

template <typename Info, typename Tag>
struct treap {
  public:

    treap() : d(1, {Info::identity(), Info::identity(), {}, 0, 0, 0, 0, false}), root(0) { }

    void insert(unsigned int idx, Info info, unsigned int priority) {
        unsigned int l, r;
        split(root, l, r, idx);
        root = merge(merge(l, new_node(info, priority)), r);
    }

    void append(Info info, unsigned int priority) {
        root = merge(root, new_node(info, priority));
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

    Info range_query(unsigned int l, unsigned int r) {
        unsigned int a, b, c;
        split(root, b, c, r + 1);
        split(b, a, b, l);
        Info ret = d[b].info;
        root = merge(merge(a, b), c);
        return ret;
    }

    void range_update(unsigned int l, unsigned int r, Tag tag) {
        unsigned int a, b, c;
        split(root, b, c, r + 1);
        split(b, a, b, l);
        d[b].lazy.apply(tag);
        root = merge(merge(a, b), c);
    }

    unsigned int size() {
        return d[root].cnt;
    }

  private:
    
    struct node {
        Info val, info; Tag lazy;
        unsigned int cnt, pri, left, right;
        bool rev;
    };

    std::vector<node> d;
    std::vector<unsigned int> del;
    unsigned int root;

    unsigned int new_node(Info val, unsigned int priority) {
        unsigned int id;
        if (del.empty()) {
            id = d.size();
            d.push_back({val, val, {}, 1, priority, 0, 0, false});
        } else {
            id = del.back();
            del.pop_back();
            d[id] = {val, val, {}, 1, priority, 0, 0, false};
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
            if (d[id].left) d[d[id].left].rev ^= true;
            if (d[id].right) d[d[id].right].rev ^= true;
        }
        if (d[id].lazy == Tag()) return;
        d[id].lazy = Tag();
        d[id].val.apply(d[id].lazy, 1);
        d[id].info.apply(d[id].lazy, d[id].cnt);
        if (d[id].left) d[d[id].left].lazy.apply(d[id].lazy);
        if (d[id].right) d[d[id].right].lazy.apply(d[id].lazy);
    }

    void update_node(unsigned int id) {
        if (!id) return;
        d[id].cnt = 1 + d[d[id].left].cnt + d[d[id].right].cnt;
        d[id].info = Info::merge(d[id].val, Info::merge(d[d[id].left].info, d[d[id].right].info));
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
