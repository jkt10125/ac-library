#ifndef ATCODER_STACK_HPP
#define ATCODER_STACK_HPP 1

#include <vector>
#include <array>

namespace atcoder {

template <class S, S (*op)(S, S), S (*e)()> struct stack {

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