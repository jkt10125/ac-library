#ifndef ATCODER_CONVEX_HULL_HPP
#define ATCODER_CONVEX_HULL_HPP 1

#include "point.hpp"

#include <vector>
#include <algorithm>

namespace atcoder {

std::vector<point> convex_hull(std::vector<point> pts) {
    if (pts.size() <= 1) return pts;
    sort(pts.begin(), pts.end());
    std::vector<point> H(pts.size() + 1);
    int s = 0, t = 0;
    for (int i = 2; i--; s = --t, std::reverse(pts.begin(), pts.end())) {
        for (point &p : pts) {
            while (t > s + 1 && H[t - 2].cross(H[t - 1], p) <= 0) {
                t--;
            }
            H[t++] = p;
        }
    }
    return {H.begin(), H.begin() + t - (t == 2 && H[0] == H[1])};
};

} // namespace atcoder

#endif // ATCODER_CONVEX_HULL_HPP