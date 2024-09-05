#ifndef ATCODER_RANDOM_HPP
#define ATCODER_RANDOM_HPP 1

#include <random>
#include <chrono>

namespace atcoder {

// random with uniform_int_distribution

std::mt19937_64 rng(std::chrono::steady_clock::now().time_since_epoch().count());

unsigned long long rand_int64(unsigned long long l, unsigned long long r) {
    return std::uniform_int_distribution<unsigned long long>(l, r)(rng);
}

}  // namespace atcoder

#endif  // ATCODER_RANDOM_HPP