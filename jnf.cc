#include <iostream>
#include <pair>
#include <random>
#include <vector>

#include "matrix.h"
#include "polynomial.h"
#include "rational.h"
#include "defs.h"

using namespace std;

// 1.1
// vector<vector<Rational>> arr = {
//     {6, 5, 1, 3},
//     {-3, -3, -1, -4},
//     {-3, -4, 2, -2},
//     {3, 6, 1, 7}
// };

// vector<pair<Rational, u32>> char_poly = {
//     {3, 4}
// };

// 1.2
// vector<vector<Rational>> arr = {
//     {2, -5, 1, 3},
//     {1, 7, -1, -2},
//     {1, 4, 2, -2},
//     {1, 4, -1, 1}
// };

// vector<pair<Rational, u32>> char_poly = {
//     {3, 4}
// };

// 1.3
// vector<vector<Rational>> arr = {
//     {6, -2, 2, -1},
//     {5, -1, 3, -2},
//     {-3, 2, 1, 1},
//     {-7, 6, -4, 6}
// };

// vector<pair<Rational, u32>> char_poly = {
//     {3, 4}
// };


// 2
vector<vector<Rational>> arr = {
    {5, 3, -1, 0, 0, 0},
    {0, 3, 1, 0, 0, 0},
    {0, -2, 6, 0, 0, 0},
    {-3, -1, -3, 7, 1, 0},
    {6, 1, 8, -4, 3, 0},
    {-5, -5, -4, 2, 2, 4}
};

vector<pair<Rational, u32>> char_poly = {
    {5, 4},
    {4, 2}
};

mt19937 gen(random_device{}());

int main() {
    Matrix<Rational> A(arr);
    vector<pair<Rational, u32>> cell_sizes;

    // Calculating cells
    for (auto p : char_poly) {
        auto B = A - Matrix<Rational>(A).ones() * p.first;
        vector<u32> pow_ranks(A.n + 2);

        for (u32 size = 0; size <= A.n + 1; ++size) {
            pow_ranks[size] = rk(pow(B, size));
        }

        for (u32 size = 1; size <= p.second; ++size) {
            u32 cell_cnt = pow_ranks[size - 1] + pow_ranks[size + 1] - 2 * pow_ranks[size];
            if (cell_cnt) {
                for (u32 i = 0; i < cell_cnt; ++i) {
                    cell_sizes.push_back({ p.first, size});
                }
                cout << "lambda: " << cell_sizes.back().first
                     << " cell size: " << size
                     << " num cells: " << cell_cnt << "\n";
            }
        }
    }

    vector<Matrix<Rational>> jnf_base;

    // This part doesn't always work but I got a correct basis at some point and used it
    // Calculating basis
    for (auto p : cell_sizes) {
        auto B = A - Matrix<Rational>(A).ones() * p.first;

        auto ker = fss(pow(B, p.second));

        vector<Rational> coefficients(ker.size());

        while (true) {
            coefficients[gen() % coefficients.size()] += 1;
            Matrix<Rational> base(A.n, 1);
            vector<Matrix<Rational>> current_basis;

            for (u32 i = 0; i < ker.size(); ++i) {
                base = base + ker[i] * coefficients[i];
            }

            for (u32 i = 0; i < p.second; ++i) {
                current_basis.push_back(base);
                base = B * base;
            }

            if (independent_vectors(current_basis)) {
                for (auto v : current_basis) {
                    cout << v << "---\n";
                    jnf_base.push_back(v);
                }
                break;
            }
        }
    }

    return 0;
}