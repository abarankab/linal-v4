#ifndef __MATRIX_H
#define __MATRIX_H

#include <pair>
#include <vector>
#include <stdexcept>

#include "defs.h"

using namespace std;

template<typename T>
struct Matrix {
    u32 n;
    u32 m;
    vector<vector<T>> data;

    Matrix(u32 n, u32 m) : n(n), m(m), data(vector<vector<T>> (n, vector<T>(m, 0))) {}

    Matrix(u32 n, u32 m, T x) : n(n), m(m), data(vector<vector<T>> (n, vector<T>(m, 0))) {
        for (u32 i = 0; i < n; ++i) {
            for (u32 j = 0; j < m; ++j) {
                if (i == j) data[i][j] = x;
            }
        }
    }

    Matrix(const vector<vector<T>>& data) : n(data.size()), m(data[0].size()), data(data) {}

    vector<T>& operator[](u32 id) {
        return data[id];
    }

    vector<T> operator[](u32 id) const {
        return data[id];
    }

    Matrix<T>& zeros() {
        for (u32 i = 0; i < n; ++i) {
            for (u32 j = 0; j < m; ++j) {
                data[i][j] = 0;
            }
        }

        return (*this);
    }

    Matrix<T>& ones() {
        for (u32 i = 0; i < n; ++i) {
            for (u32 j = 0; j < m; ++j) {
                data[i][j] = (i == j) ? 1 : 0;
            }
        }

        return (*this);
    }

    T trace() {
        if (n != m) {
            throw invalid_argument("Trace of a non-square matrix is undefined");
        }

        T R = 0;
        for (u32 i = 0; i < n; ++i) {
            R += data[i][i];
        }

        return R;
    }

    Matrix<T>& transpose() {
        vector<vector<T>> new_data(m, vector<T>(n));
        for (u32 i = 0; i < n; ++i) {
            for (u32 j = 0; j < m; ++j) {
                new_data[j][i] = data[i][j];
            }
        }
        swap(n, m);
        data = new_data;
        return (*this);
    }

    Matrix<T>& transposed() {
        Matrix<T> R = (*this);
        return R.transpose();
    }
};

template<typename T>
Matrix<T> operator*(const Matrix<T>& A, const Matrix<T>& B) {
    if (A.m != B.n) {
        throw invalid_argument("Matrix multiplication incorrect dimensions");
    }

    Matrix<T> R(A.n, B.m);
    
    for (u32 i = 0; i < R.n; ++i) {
        for (u32 j = 0; j < R.m; ++j) {
            for (u32 k = 0; k < A.m; ++k) {
                R[i][j] += A[i][k] * B[k][j];
            }
        }
    }

    return R;
}

template<typename T>
Matrix<T>& operator*=(Matrix<T>& A, const Matrix<T>& B) {
    return A = (A * B);
}

template<typename T>
Matrix<T> operator+(const Matrix<T>& A, const Matrix<T>& B) {
    if (A.n != B.n || A.m != B.m) {
        throw invalid_argument("Matrix addition incorrect dimensions");
    }

    Matrix<T> R(A.n, A.m);

    for (u32 i = 0; i < R.n; ++i) {
        for (u32 j = 0; j < R.m; ++j) {
            R[i][j] = A[i][j] + B[i][j];
        }
    }

    return R;
}

template<typename T>
Matrix<T>& operator+=(Matrix<T>& A, const Matrix<T>& B) {
    return A = (A + B);
}

template<typename T>
Matrix<T> operator-(const Matrix<T>& A, const Matrix<T>& B) {
    if (A.n != B.n || A.m != B.m) {
        throw invalid_argument("Matrix addition incorrect dimensions");
    }

    Matrix<T> R(A.n, A.m);

    for (u32 i = 0; i < R.n; ++i) {
        for (u32 j = 0; j < R.m; ++j) {
            R[i][j] = A[i][j] - B[i][j];
        }
    }

    return R;
}

template<typename T>
Matrix<T>& operator-=(Matrix<T>& A, const Matrix<T>& B) {
    return A = (A - B);
}

template<typename T>
Matrix<T> operator*(const Matrix<T>& A, const T& B) {
    Matrix<T> R(A);

    for (u32 i = 0; i < R.n; ++i) {
        for (u32 j = 0; j < R.m; ++j) {
            R[i][j] *= B;
        }
    }

    return R;
}

template<typename T>
Matrix<T>& operator*=(Matrix<T>& A, const T& B) {
    return A = (A * B);
}

template <typename T>
ostream& operator<<(ostream& os, Matrix<T> A) {
    for (u32 i = 0; i < A.n; ++i) {
        for (u32 j = 0; j < A.m; ++j) {
            os << A[i][j] << " ";
        }
        os << "\n";
    }
    return os;
}

template<typename T>
Matrix<T> pow(const Matrix<T>& A, u32 n) {
    if (n == 0) {
        return Matrix<T>(A).ones();
    }

    Matrix<T> R = A;
    for (u32 i = 0; i < n - 1; ++i) {
        R *= A;
    }

    return R;
}

template<typename T>
Matrix<T> row_reduce(const Matrix<T>& A) {
    Matrix<T> R = A;
    for (u32 i = 0; i < R.n; ++i) {
        u32 fpos = R.m - 1;
        u32 fid = i;

        for (u32 j = i; j < R.n; ++j) {
            for (u32 k = 0; k < R.m; ++k) {
                if (R[j][k] != 0 && fpos > k) {
                    fpos = k;
                    fid = j;
                    break;
                }
            }
        }

        if (fpos != R.m - 1) {
            swap(R.data[i], R.data[fid]);
        }

        u32 pos = 0;
        bool found = false;

        for (u32 k = 0; k < R.m; ++k) {
            if (R[i][k] != 0) {
                found = true;
                pos = k;
                break;
            }
        }

        if (!found) break;

        T multipler = R[i][pos];
        for (u32 k = pos; k < R.m; ++k) {
            R[i][k] /= multipler;
        }

        for (u32 j = i + 1; j < R.n; ++j) {
            T c = R[j][pos];
            for (u32 k = 0; k < R.m; ++k) {
                R[j][k] -= R[i][k] * c;
            }
        }
    }

    return R;
}

template<typename T>
u32 rk(const Matrix<T>& A) {
    auto R = row_reduce(A);
    u32 res = 0;
    for (u32 i = 0; i < R.n; ++i) {
        bool is_null = true;
        for (u32 j = 0; j < R.m; ++j) {
            is_null &= (R[i][j] == 0);
        }
        if (!is_null) ++res;
    }
    return res;
}

template<typename T>
vector<Matrix<T>> fss(const Matrix<T>& A) {
    auto R = row_reduce(A);
    vector<u32> free_positions;
    vector<u32> main_positions;
    for (u32 i = 0, offset = 0; i < R.n; ++i) {
        while (i + offset != R.m && i != R.n && R[i][i + offset] != 1) {
            free_positions.push_back(i + offset++);
        }
        if (i + offset != R.m && i != R.n) {
            main_positions.push_back(i + offset++);
        } else {
            break;
        }
    }

    vector<Matrix<T>> result(free_positions.size(), Matrix<T>(R.m, 1));
    for (u32 i = 0; i < free_positions.size(); ++i) {
        for (u32 j : main_positions) {
            if (j > free_positions[i]) break;

            result[i][j][0] = R[j][free_positions[i]] * -1;
        }
        result[i][free_positions[i]][0] = 1;
    }

    return result;
}

template<typename T>
bool independent_vectors(vector<Matrix<T>> vectors) {
    Matrix<T> A(vectors.size(), vectors[0].n);
    for (u32 i = 0; i < vectors.size(); ++i) {
        for (u32 k = 0; k < vectors[0].n; ++k) {
            A[i][k] = vectors[i][k][0];
        }
    }

    return rk(A) == A.n;
}

#endif
