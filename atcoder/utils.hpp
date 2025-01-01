#ifndef ATCODER_UTILS_HPP
#define ATCODER_UTILS_HPP 1

#include <vector>
#include <random>
#include <algorithm>
#include <chrono>

namespace atcoder {

// Return a random permutation of size n. (0-indexed)
std::vector<int> rand_perm(size_t n) {
    std::mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count());
    std::vector<int> p(n);
    std::iota(p.begin(), p.end(), 0);
    std::shuffle(p.begin(), p.end(), rng);
    return p;
}

}

#endif  // ATCODER_UTILS_HPP