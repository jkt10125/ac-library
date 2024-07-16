#ifndef ATCODER_POINT_HPP
#define ATCODER_POINT_HPP 1

namespace atcoder {

using S = int;
// template <class T>
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

} // namespace atcoder

#endif // ATCODER_POINT_HPP