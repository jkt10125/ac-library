#ifndef ATCODER_STACK_HPP
#define ATCODER_STACK_HPP 1

#include <vector>
#include <array>

namespace atcoder {

#if __cplusplus >= 201703L

template <class S, auto op, auto e> struct stack {
    static_assert(std::is_convertible_v<decltype(op), std::function<S(S, S)>>,
                  "op must work as S(S, S)");
    static_assert(std::is_convertible_v<decltype(e), std::function<S()>>,
                  "e must work as S()");

#else

template <class S, S (*op)(S, S), S (*e)()> struct stack {

#endif

  public:
    void push(S x) {
        d.push_back({x, op(x, prod())});
    }

    void pop() {
        d.pop_back();
    }

    S top() {
        return d.empty() ? e() : d.back()[0];
    }

    S prod() {
        return d.empty() ? e() : d.back()[1];
    }

    int size() {
        return int(d.size());
    }

    bool empty() {
        return d.empty();
    }

    void clear() {
        d.clear();
    }

  private:
    std::vector<std::array<S, 2>> d;
};


}  // namespace atcoder

#endif  // ATCODER_STACK_HPP