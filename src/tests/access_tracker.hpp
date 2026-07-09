#pragma once

template<typename T>
struct CountedValue {
    T value{};
    mutable size_t access_count = 0;

    CountedValue() = default;
    CountedValue(T v) : value(v) {}  // Add this
    operator T() const {
        ++access_count;
        return value;
    }

    CountedValue& operator=(const T& v) {
        value = v;
        return *this;
    }
};