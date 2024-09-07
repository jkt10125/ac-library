#ifndef ATCODER_MISC_ALGO_HPP
#define ATCODER_MISC_ALGO_HPP 1

#include <vector>
#include <array>

namespace atcoder {

std::vector<int> articulation_points(const std::vector<std::vector<int>>& A) {
    int n = A.size();
    std::vector<int> ord(n, -1), low(n), ret;
    int now_ord = 0;
    auto dfs = [&](auto self, int v, int p) -> void {
        ord[v] = low[v] = now_ord++;
        int cnt = 0;
        bool is_articulation = false;
        for (int to : A[v]) {
            if (to == p) continue;
            if (ord[to] == -1) {
                cnt++;
                self(self, to, v);
                low[v] = std::min(low[v], low[to]);
                if (p != -1 && ord[v] <= low[to]) {
                    is_articulation = true;
                }
            } else {
                low[v] = std::min(low[v], ord[to]);
            }
        }
        if (p == -1 && cnt > 1) is_articulation = true;
        if (is_articulation) ret.push_back(v);
    };
    for (int i = 0; i < n; i++) {
        if (ord[i] == -1) dfs(dfs, i, -1);
    }
    return ret;
}

std::vector<std::array<int, 2>> bridges(const std::vector<std::vector<int>>& A) {
    int n = A.size();
    std::vector<int> ord(n, -1), low(n);
    std::vector<std::array<int, 2>> ret;
    int now_ord = 0;
    auto dfs = [&](auto self, int v, int p) -> void {
        ord[v] = low[v] = now_ord++;
        for (int to : A[v]) {
            if (to == p) continue;
            if (ord[to] == -1) {
                self(self, to, v);
                low[v] = std::min(low[v], low[to]);
                if (ord[v] < low[to]) {
                    ret.push_back({v, to});
                }
            } else {
                low[v] = std::min(low[v], ord[to]);
            }
        }
    };
    for (int i = 0; i < n; i++) {
        if (ord[i] == -1) dfs(dfs, i, -1);
    }
    return ret;
}

std::vector<int> topo_sort(const std::vector<std::vector<int>>& A) {
    int n = A.size();
    std::vector<int> ret, in(n);
    for (int i = 0; i < n; i++) {
        for (int j : A[i]) {
            in[j]++;
        }
    }
    std::vector<int> que;
    for (int i = 0; i < n; i++) {
        if (in[i] == 0) {
            que.push_back(i);
        }
    }
    while (!que.empty()) {
        int v = que.back();
        que.pop_back();
        ret.push_back(v);
        for (int u : A[v]) {
            in[u]--;
            if (in[u] == 0) {
                que.push_back(u);
            }
        }
    }
    return ret;
}

long long tsp(const std::vector<std::vector<long long>>& G) {
    unsigned int n = G.size();
    std::vector<std::vector<long long>> dp(1 << n, std::vector<long long>(n, 1LL << 60));
    dp[1][0] = 0;
    for (int S = 1; S < (1 << n); S += 2) {
        for (int v = 0; v < n; v++) {
            if (!(S >> v & 1)) continue;
            for (int u = 0; u < n; u++) {
                if (S >> u & 1) {
                    dp[S][v] = std::min(dp[S][v], dp[S ^ (1 << v)][u] + G[u][v]);
                }
            }
        }
    }
    return dp[(1 << n) - 1][0];
}

}  // namespace atcoder

#endif  // ATCODER_MISC_ALGO_HPP