#ifndef ATCODER_QUEUE_HPP
#define ATCODER_QUEUE_HPP 1

#include "stack.hpp"

#include <array>

namespace atcoder {

template <class S, S (*op)(S, S), S (*e)()> struct queue {

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