#ifndef ATCODER_TREE_ANCESTOR_HPP
#define ATCODER_TREE_ANCESTOR_HPP 1

#include <vector>
#include <algorithm>

namespace atcoder {

struct treeancestor {
  public:

    treeancestor() { }
    treeancestor(int n) : _n(n), time(0) {
        sz.resize(n);       // size of the subtree
        in.resize(n);       // in-time of dfs
        out.resize(n);      // out-time of dfs
        seq.resize(n);      // the dfs order
        lvl.resize(n);      // distance from root
        par.resize(n);      // parent of node
        top.resize(n);      // top node of HLD-thread
        adj.assign(n, {});  // adjacency list
    }

    void add_edge(int u, int v) {
        adj[u].push_back(v);
        adj[v].push_back(u);
    }

    void apply(int root) {
        top[root] = root;
        dfs1(root, -1), dfs2(root);
    }

    int lca(int x, int y) {
        while (top[x] != top[y]) {
            if (lvl[top[x]] > lvl[top[y]]) { x = par[top[x]]; }
            else { y = par[top[y]]; }
        }
        return (lvl[x] < lvl[y]) ? x : y;
    }

    int dist(int x, int y) {
        return lvl[x] + lvl[y] - 2 * lvl[lca(x, y)];
    }

    int jump(int x, int k) {
        if (lvl[x] < k) { return -1; }
        int depth = lvl[x] - k;
        while (lvl[top[x]] > depth) {
            x = par[top[x]];
        }
        return seq[in[x] - lvl[x] + depth];
    }
    
    // x ancestor of y ?
    bool is_ancestor(int x, int y) {
        return in[x] < in[y] && in[y] < out[x];
    }
    
    // tree rooted at r
    int rooted_parent(int r, int x) {
        if (r == x) { return x; }
        if (!is_ancestor(x, r)) {
            return par[x];
        }
        auto itr = std::upper_bound(adj[r].begin(), adj[r].end(), x, [&] (int a, int b) { return in[a] < in[b]; });
        return *(std::prev(itr));
    }

    int rooted_size(int r, int x) {
        if (r == x) { return _n; }
        if (!is_ancestor(x, r)) {
            return sz[x];
        }
        return _n - sz[rooted_parent(r, x)];
    }

    int rooted_lca(int a, int b, int c) {
        return lca(a, b) ^ lca(b, c) ^ lca(c, a);
    }

  private:
    int _n, time;
    std::vector<int> sz, top, lvl, par, in, out, seq;
    std::vector<std::vector<int>> adj;

    void dfs1(int x, int p) {
        if ((par[x] = p) == -1) { lvl[x] = 0; }
        else {
            adj[x].erase(std::find(adj[x].begin(), adj[x].end(), p));
            lvl[x] = lvl[p] + 1;
        }
        sz[x] = 1;
        for (int &it : adj[x]) {
            dfs1(it, x);
            sz[x] += sz[it];
            if (sz[it] > sz[adj[x][0]]) {
                std::swap(it, adj[x][0]);
            }
        }
    }

    void dfs2(int x) {
        in[x] = time++, seq[in[x]] = x;
        for (int it : adj[x]) {
            top[it] = (it == adj[x].front()) ? top[x] : it;
            dfs2(it);
        }
        out[x] = time;
    }
};

} // namespace atcoder

#endif  // ATCODER_TREE_ANCESTOR_HPP
