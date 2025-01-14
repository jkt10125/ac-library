#ifndef ATCODER_TREE_ANCESTOR_OP_HPP
#define ATCODER_TREE_ANCESTOR_OP_HPP 1

#include <vector>
#include <array>
#include <functional>
#include <cassert>

namespace atcoder {

template <class S, S (*op)(S, S), S (*e)()> struct treeancestor_op {

  public:
    
    treeancestor_op(const std::vector<std::vector<std::pair<int, S>>>& A, int root) : _n(A.size()), log(0) {
        while ((1 << log) < _n) log++;
        lvl.resize(_n);
        lift.resize(log, std::vector<int>(_n));
        d.resize(log, std::vector<S>(_n));

        dfs(root, -1, A);
    }

    treeancestor_op(const std::vector<int>& par,
                    const std::vector<S>& W) {
        _n = par.size();
        int root = 0;
        while (root < _n && par[root] != -1) { root++; }
        
        assert(root < _n);
        for (int i = root + 1; i < _n; i++) 
            assert(par[i] != -1);

        std::vector<std::vector<std::pair<int, S>>> A(_n);
        for (int i = 0; i < _n; i++) {
            if (i == root) continue;
            A[par[i]].push_back({i, W[i]});
            A[i].push_back({par[i], W[i]});
        }
        *this = treeancestor_op(A, root);
    }

    int kth_ancestor(int x, int k) {
        for (int i = log - 1; i >= 0; i--) {
            if (k & (1 << i)) { x = lift[i][x]; }
        }
        return x;
    }

    int lca(int x, int y) {
        if (lvl[x] < lvl[y]) { std::swap(x, y); }
        x = kth_ancestor(x, lvl[x] - lvl[y]);
        if (x == y) { return x; }
        for (int i = log - 1; i >= 0; i--) {
            if (lift[i][x] != lift[i][y]) {
                x = lift[i][x];
                y = lift[i][y];
            }
        }
        return lift[0][x];
    }

    S query(int x, int y) {
        int l = lca(x, y);
        S res = e();
        for (int i = log - 1; i >= 0; i--) {
            if (lvl[x] - (1 << i) >= lvl[l]) {
                res = op(res, d[i][x]);
                x = lift[i][x];
            }
            if (lvl[y] - (1 << i) >= lvl[l]) {
                res = op(res, d[i][y]);
                y = lift[i][y];
            }
        }
        return res;
    }

  private:
    int _n, log;
    std::vector<int> lvl;
    std::vector<std::vector<S>> d;
    std::vector<std::vector<int>> lift;

    void dfs(int x, int p, const std::vector<std::vector<std::pair<int, S>>>& A) {
        lvl[x] = (p != -1) ? lvl[p] + 1 : 0;
        lift[0][x] = p;
        for (int i = 1; i < log; i++) {
            lift[i][x] = (lift[i - 1][x] != -1) ? lift[i - 1][lift[i - 1][x]] : -1;
            d[i][x] = op(d[i - 1][x], (lift[i - 1][x] != -1) ? d[i - 1][lift[i - 1][x]] : e());
        }
        for (auto [y, w] : A[x]) {
            if (y != p) {
                d[0][y] = w;
                dfs(y, x, A);
            }
        }
    }
};

} // namespace atcoder

#endif  // ATCODER_TREE_ANCESTOR_OP_HPP
