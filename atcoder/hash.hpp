#ifndef ATCODER_HASH_HPP
#define ATCODER_HASH_HPP 1

#include <vector>
#include <string>
#include <cstring>

namespace atcoder {

// Reference: https://richardstartin.github.io/posts/xxhash
struct XXHash64 {
  public:
    XXHash64() : prime1(11400714785074694791ULL),
                 prime2(14029467366897019727ULL),
                 prime3(1609587929392839161ULL),
                 prime4(9650029242287828579ULL),
                 prime5(2870177450012600261ULL) { }
    
    uint64_t operator()(const std::vector<uint8_t>& data, uint64_t seed) {
        int len = data.size();
        auto it = data.begin();
        auto end = data.end();
        uint64_t hash64;

        if (len >= 32) {
            auto limit = it + (len - 32);
            uint64_t v1 = seed + prime1 + prime2;
            uint64_t v2 = seed + prime2;
            uint64_t v3 = seed;
            uint64_t v4 = seed - prime1;

            while (it <= limit) {
                v1 += get_uint64(it) * prime2;
                v1 = rotate_left(v1, 31);
                v1 *= prime1;

                v2 += get_uint64(it + 8) * prime2;
                v2 = rotate_left(v2, 31);
                v2 *= prime1;

                v3 += get_uint64(it + 16) * prime2;
                v3 = rotate_left(v3, 31);
                v3 *= prime1;

                v4 += get_uint64(it + 24) * prime2;
                v4 = rotate_left(v4, 31);
                v4 *= prime1;

                it += 32;
            }

            hash64 = rotate_left(v1, 1) + rotate_left(v2, 7) + rotate_left(v3, 12) + rotate_left(v4, 18);
        } else {
            hash64 = seed + prime5;
        }

        hash64 += len;

        while (std::distance(it, end) >= 8) {
            hash64 += get_uint64(it) * prime3;
            hash64 = rotate_left(hash64, 27) * prime4;
            std::advance(it, 8);
        }

        while (std::distance(it, end) >= 4) {
            hash64 += get_uint32(it) * prime3;
            hash64 = rotate_left(hash64, 23) * prime4;
            std::advance(it, 4);
        }

        while (it != end) {
            hash64 += (*it) * prime5;
            hash64 = rotate_left(hash64, 11) * prime1;
            ++it;
        }

        hash64 ^= hash64 >> 33;
        hash64 *= prime2;
        hash64 ^= hash64 >> 29;
        hash64 *= prime3;
        hash64 ^= hash64 >> 32;

        return hash64;
    }

  private:
    uint64_t prime1, prime2, prime3, prime4, prime5;

    uint64_t get_uint64(std::vector<uint8_t>::const_iterator it) const {
        uint64_t result;
        std::memcpy(&result, &(*it), sizeof(result));
        return result;
    }

    uint32_t get_uint32(std::vector<uint8_t>::const_iterator it) const {
        uint32_t result;
        std::memcpy(&result, &(*it), sizeof(result));
        return result;
    }

    uint64_t rotate_left(uint64_t x, int r) const {
        return (x << r) | (x >> (64 - r));
    }
};

};  // namespace atcoder

#endif  // ATCODER_HASH_HPP