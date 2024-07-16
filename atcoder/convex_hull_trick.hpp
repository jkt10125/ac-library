#ifndef ATCODER_CONVEX_HULL_TRICK_HPP
#define ATCODER_CONVEX_HULL_TRICK_HPP 1

#include <set>
#include <array>
#include <vector>
#include <limits>
#include <iterator>

namespace atcoder {

namespace internal {

template <class S>
struct line {
    mutable S m, c, p;
    bool isline;
    bool operator < (const line& l) const {
        return (l.isline) ? (m < l.m) : (p < l.p);
    }
};

} // namespace internal

using S = long long;
struct convex_hull_trick : std::multiset<internal::line<S>> {
  public:
	convex_hull_trick() : mode(1) { } // default mode is max

	inline void flip_mode() { // flip mode between max and min
        mode *= -1;
    } 
	
    void add (S slope, S intercept) {
    	auto z = insert({mode * slope, mode * intercept, 0, 1});
        auto y = z++;
        auto x = y;
    	while (isect(y, z)) { z = erase(z); }
    	if (x != begin() && isect(--x, y)) {
            isect(x, y = erase(y));
        }
    	while ((y = x) != begin() && (--x)->p >= y->p) {
            isect(x, erase(y));
        }
	}
	
    S query (S x) {
    	// assert(!empty());
		auto l = *lower_bound({0, 0, x, 0});
    	return mode * (l.m * x + l.c);
	}
  
  private:
    int mode;

	S div (S a, S b) {
        return (a / b - ((a ^ b) < 0 && a % b));
    }
	
	bool isect (iterator x, iterator y){
    	if (y == end()) {
            x->p = std::numeric_limits<S>::max();
            return false;
        }
    	x->p = (x->m == y->m) ? 
                    (x->c > y->c) ? 
                            std::numeric_limits<S>::max() : 
                            std::numeric_limits<S>::min() : 
                    div ((y->c - x->c), (x->m - y->m));
        return (x->p >= y->p);
	}
};

} // namespace atcoder

#endif // ATCODER_CONVEX_HULL_TRICK_HPP
