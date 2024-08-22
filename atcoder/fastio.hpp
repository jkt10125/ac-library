#ifndef ATCODER_FASTIO_HPP
#define ATCODER_FASTIO_HPP 1

#include "internal_type_traits.hpp"

#include <string>
#include <vector>
#include <algorithm>

namespace atcoder {

void fast_scan(std::string& x) {
    x.clear();
    char c = getchar_unlocked();
    while (isspace(c)) {
        c = getchar_unlocked();
    }
    for (; !isspace(c); c = getchar_unlocked()) {
        x += c;
    }
}

template <class T, std::enable_if_t<internal::is_truely_char<T>::value, std::nullptr_t> = nullptr>
void fast_scan(T& x) {
    char c = getchar_unlocked();
    while (isspace(c)) {
        c = getchar_unlocked();
    }
    x = c;
}

void fast_scan(float& x) {
    std::string s;
    fast_scan(s);
    x = std::stof(std::move(s));
}

void fast_scan(double& x) {
    std::string s;
    fast_scan(s);
    x = std::stod(std::move(s));
}

void fast_scan(long double& x) {
    std::string s;
    fast_scan(s);
    x = std::stold(std::move(s));
}

template <class T, std::enable_if_t<internal::is_truely_integral<T>::value, std::nullptr_t> = nullptr>
void fast_scan(T& x) {
    using uT = internal::to_unsigned_t<T>;

    uT y = 0;
    bool neg = false;
    char c = getchar_unlocked();
    while (isspace(c)) {
        c = getchar_unlocked();
    }
    if (c == '-') {
        neg = true;
        c = getchar_unlocked();
    }
    for (; !isspace(c); c = getchar_unlocked()) {
        y = y * 10 + (c - '0');
    }
    x = neg ? -y : y;
}

template <class T, class U>
void fast_scan(std::pair<T, U>& x); // forward declaration for the case of std::vector<std::pair ...>

template <class T, std::enable_if_t<internal::is_iterable<T>::value, std::nullptr_t> = nullptr>
void fast_scan(T& x) {
    for (auto& y : x) {
        fast_scan(y);
    }
}

template <class T, class U>
void fast_scan(std::pair<T, U>& x) {
    fast_scan(x.first);
    fast_scan(x.second);
}

// ------------------------------------------------------------------------------------------

void fast_print(const std::string& x) {
    for (char c : x) {
        putchar_unlocked(c);
    }
}

template <class T, std::enable_if_t<internal::is_truely_char<T>::value, std::nullptr_t> = nullptr>
void fast_print(const T& x) {
    putchar_unlocked(x);
}

template <class T, std::enable_if_t<std::is_floating_point<T>::value, std::nullptr_t> = nullptr>
void fast_print(const T& x) {
    std::string s = std::to_string(x);
    fast_print(s);
}

template <class T, std::enable_if_t<internal::is_truely_integral<T>::value, std::nullptr_t> = nullptr>
void fast_print(const T& x) {
    using uT = internal::to_unsigned_t<T>;

    if (x == 0) {
        putchar_unlocked('0');
    }
    else {
        uT y = x < 0 ? -x : x;
        if (x < 0) {
            putchar_unlocked('-');
        }
        std::string s;
        while (y > 0) {
            s += char('0' + (y % 10));
            y /= 10;
        }
        std::reverse(s.begin(), s.end());
        fast_print(s);
    }
}

template <class T, class U>
void fast_print(const std::pair<T, U>& x, char separator = ' '); // forward declaration for the case of std::vector<std::pair ...>

template <class T, std::enable_if_t<internal::is_iterable<T>::value, std::nullptr_t> = nullptr>
void fast_print(const T& x, char separator = ' ') {
    for (auto it = x.cbegin(); it != x.cend(); ++it) {
        if (it != x.cbegin()) {
            putchar_unlocked(separator);
        }
        fast_print(*it);
    }
}

template <class T, class U>
void fast_print(const std::pair<T, U>& x, char separator) {
    fast_print(x.first);
    putchar_unlocked(separator);
    fast_print(x.second);
}

}  // namespace atcoder

#endif      // ATCODER_FASTIO_HPP