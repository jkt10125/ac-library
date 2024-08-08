#ifndef ATCODER_QUEUE_HPP
#define ATCODER_QUEUE_HPP 1

#include "stack.hpp"

#include <array>

namespace atcoder {

#if __cplusplus >= 201703L

template <class S, auto op, auto e> struct queue {
    static_assert(std::is_convertible_v<decltype(op), std::function<S(S, S)>>,
                  "op must work as S(S, S)");
    static_assert(std::is_convertible_v<decltype(e), std::function<S()>>,
                  "e must work as S()");

#else

template <class S, S (*op)(S, S), S (*e)()> struct queue {

#endif

  public:
    
    void push(S x) {
        entry.push(x);
    }

    void pop() {
        stack_transfer();
        exit.pop();
    }

    bool empty() {
        return entry.empty() && exit.empty();
    }

    int size() {
        return int(entry.size() + exit.size());
    }

    S front() {
        stack_transfer();
        return exit.top();
    }

    S prod() {
        return op(entry.prod(), exit.prod());
    }

  private:
    stack<S, op, e> entry, exit;
    
    void stack_transfer() {
        if (exit.empty()) {
            while (!entry.empty()) {
                exit.push(entry.top());
                entry.pop();
            }
        }
    }

};

} // namespace atcoder

#endif // ATCODER_QUEUE_HPP