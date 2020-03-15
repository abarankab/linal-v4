#ifndef __RATIONAL_H
#define __RATIONAL_H

#include <iostream>
#include <numeric>

#include "defs.h"

using namespace std;

struct Rational {
    i64 p;
    i64 q;

    Rational(i64 x = 0) : p(x), q(1) {}

    Rational(i64 p, i64 q) : p(p), q(q) {}

    Rational& compact() {
        i64 g = gcd(p, q);
        p /= g;
        q /= g;

        if (q < 0) {
            p *= -1;
            q *= -1;
        }

        if (p == 0) q = 1;

        return (*this);
    }
};

bool operator==(const Rational A, const Rational B) {
    return (A.p == B.p) && (A.q == B.q);
}

bool operator!=(const Rational A, const Rational B) {
    return !(A == B);
}

bool operator<(const Rational A, const Rational B) {
    return (A.p * B.q) < (B.p * A.q);
}

bool operator<=(const Rational A, const Rational B) {
    return (A < B) || (A == B);
}

Rational operator*(const Rational A, const Rational B) {
    return Rational(A.p * B.p, A.q * B.q).compact();
}

Rational& operator*=(Rational& A, const Rational B) {
    return A = (A * B);
}

Rational operator+(const Rational A, const Rational B) {
    return Rational(A.p * B.q + B.p * A.q, A.q * B.q).compact();
}

Rational& operator+=(Rational& A, const Rational B) {
    return A = (A + B);
}

Rational operator-(const Rational A, const Rational B) {
    return Rational(A.p * B.q - B.p * A.q, A.q * B.q).compact();
}

Rational& operator-=(Rational& A, const Rational B) {
    return A = (A - B);
}

Rational operator/(const Rational A, const Rational B) {
    return Rational(A.p * B.q, A.q * B.p).compact();
}

Rational& operator/=(Rational& A, const Rational B) {
    return A = (A / B);
}

ostream& operator<<(ostream& os, const Rational A) {
    if (A.q == 1) {
        os << A.p;
    } else {
        os << A.p << "/" << A.q;
    }
    return os;
}

istream& operator>>(istream& is, Rational& A) {
    i32 x;
    is >> x;
    A = x;
    return is;
}

#endif
