#ifndef ATCODER_RANDOM_HPP
#define ATCODER_RANDOM_HPP 1

#include <random>
#include <chrono>
#include <vector>
#include <algorithm>

namespace atcoder {

// random with uniform_int_distribution

std::mt19937_64 rng64(std::chrono::steady_clock::now().time_since_epoch().count());
std::mt19937 rng32(std::chrono::steady_clock::now().time_since_epoch().count());

unsigned long long rand_int64(unsigned long long l, unsigned long long r) {
    return std::uniform_int_distribution<unsigned long long>(l, r)(rng64);
}

unsigned rand_int32(unsigned l, unsigned r) {
    return std::uniform_int_distribution<unsigned>(l, r)(rng32);
}

// Return a random permutation of size n. (0-indexed)
std::vector<int> rand_perm(size_t n) {
    std::vector<int> p(n);
    std::iota(p.begin(), p.end(), 0);
    std::shuffle(p.begin(), p.end(), rng32);
    return p;
}

}  // namespace atcoder

#endif  // ATCODER_RANDOM_HPP