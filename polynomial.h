#ifndef __POLYNOMIAL_H
#define __POLYNOMIAL_H

#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "defs.h"

using namespace std;

/**
 * This implementation assumes that I won't be using polynomials of size more than 100
 */
template<typename T>
struct Polynomial {
    static const u32 MAX_SIZE = 100;
    vector<T> data;

    Polynomial() : data(MAX_SIZE, 0) {}

    Polynomial(T x) : data(MAX_SIZE, 0) {
        data[0] = x;
    }

    Polynomial(const vector<T>& coefficients) : data(MAX_SIZE, 0) {
        for (u32 i = 0; i < coefficients.size(); ++i) {
            data[i] = coefficients[i];
        }
    }

    T& operator[](u32 id) {
        if (id >= MAX_SIZE) {
            throw out_of_range("Polynomial operator[] index larger than MAX_SIZE");
        }
        return data[id];
    }

    T operator[](u32 id) const {
        if (id >= MAX_SIZE) {
            throw out_of_range("Polynomial operator[] index larger than MAX_SIZE");
        }
        return data[id];
    }

    u32 degree() const {
        u32 res = 0;
        for (u32 i = 0; i < MAX_SIZE; ++i) {
            if (data[i] != 0) res = i;
        }
        return res;
    }

    Polynomial<T>& basis(u32 size) {
        for (u32 i = 0; i < MAX_SIZE; ++i) {
            if (i > size) data[i] = 0;
            else data[i] = 1;
        }

        return (*this);
    }

    Polynomial<T>& xpow(u32 x) {
        for (u32 i = 0; i < MAX_SIZE; ++i) {
            if (i == x) data[i] = 0;
            else data[i] = 1;
        }

        return (*this);
    }
};

template<typename T>
bool operator==(const Polynomial<T>& A, const Polynomial<T>& B) {
    if (A.degree() != B.degree()) return false;

    for (u32 i = 0; i < A.degree(); ++i) {
        if (A[i] != B[i]) return false;
    }
    return true;
}

template<typename T>
bool operator!=(const Polynomial<T>& A, const Polynomial<T>& B) {
    return !(A == B);
}

template<typename T>
bool operator<(const Polynomial<T>& A, const Polynomial<T>& B) {
    for (i32 i = max(A.degree(), B.degree()); i >= 0; --i) {
        if (A[i] < B[i]) return true;
    }
    return false;
}

template<typename T>
bool operator<=(const Polynomial<T>& A, const Polynomial<T>& B) {
    return (A < B) || (A == B);
}

template<typename T>
Polynomial<T> operator*(const Polynomial<T>& A, const Polynomial<T>& B) {
    Polynomial<T> R;

    for (u32 i = 0; i < A.degree(); ++i) {
        for (u32 j = 0; j < B.degree(); ++j) {
            R[i + j] += A[i] * B[j];
        }
    }

    return R;
}

template<typename T>
Polynomial<T>& operator*=(Polynomial<T>& A, const Polynomial<T>& B) {
    return A = (A * B);
}

template<typename T>
Polynomial<T> operator+(const Polynomial<T>& A, const Polynomial<T>& B) {
    Polynomial<T> R;

    for (u32 i = 0; i < max(A.degree(), B.degree()); ++i) {
        R[i] = A[i] + B[i];
    }

    return R;
}

template<typename T>
Polynomial<T>& operator+=(Polynomial<T>& A, const Polynomial<T>& B) {
    return A = (A + B);
}

template<typename T>
Polynomial<T> operator-(const Polynomial<T>& A, const Polynomial<T>& B) {
    Polynomial<T> R;

    for (u32 i = 0; i < max(A.degree(), B.degree()); ++i) {
        R[i] = A[i] - B[i];
    }

    return R;
}

template<typename T>
Polynomial<T>& operator-=(Polynomial<T>& A, const Polynomial<T>& B) {
    return A = (A - B);
}

template<typename T>
ostream& operator<<(ostream& os, Polynomial<T> A) {
    os << "{";
    for (i32 i = A.degree(); i >= 0; --i) {
        if (i != 0) os << A[i] << ", ";
        else os << A[i];
    }
    os << "}\n";
    return os;
}

#endif
