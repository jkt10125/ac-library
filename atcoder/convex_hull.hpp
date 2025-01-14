#ifndef ATCODER_CONVEX_HULL_HPP
#define ATCODER_CONVEX_HULL_HPP 1

#include <vector>
#include <algorithm>

namespace atcoder {

namespace internal {

using S = int;
struct point {
    S x, y;
    point() = default;
    point(S __x, S __y) : x(__x), y(__y) {}
    point operator - (const point &p) const { return point(x - p.x, y - p.y); }
    point operator + (const point &p) const { return point(x + p.x, y + p.y); }
    long long cross(const point &p) const { return x * 1ll * p.y - y * 1ll * p.x; }
    long long cross(const point &a, const point &b) const { return (a - *this).cross(b - *this); }
    bool operator == (const point &p) const { return (x == p.x) && (y == p.y); }
    bool operator < (const point &p) const { return (x == p.x) ? (y < p.y) : (x < p.x); }
};

} // namespace internal

std::vector<internal::point> convex_hull(std::vector<internal::point> pts) {
    if (pts.size() <= 1) return pts;
    sort(pts.begin(), pts.end());
    std::vector<internal::point> H(pts.size() + 1);
    int s = 0, t = 0;
    for (int i = 2; i--; s = --t, std::reverse(pts.begin(), pts.end())) {
        for (internal::point &p : pts) {
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