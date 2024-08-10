#include "fenwicktree.hpp"

#include <iostream>

signed main() {
    atcoder::fenwick_tree<long long> ft(10);
    ft.range_add(4, 4, 6);
    ft.range_add(6, 7, -5);

    for (int i = 0; i < 10; i++) {
        std::cout << ft.range_sum(i, i) << ' ';
    }
    std::cout << std::endl;

    std::cout << ft.range_sum(4, 7) << std::endl;

}