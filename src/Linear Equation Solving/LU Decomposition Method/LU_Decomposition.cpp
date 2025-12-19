#include <bits/stdc++.h>
using namespace std;
const double eps = 1e-9;

void printMatrix(ofstream &fout, vector<vector<double>> &a, int n, int m) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            fout << setw(10) << fixed << setprecision(3) << a[i][j] << " ";
        }
        fout << endl;
    }
    fout << endl;
}

int main() {
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    fout << "LU DECOMPOSITION (Doolittle)\n\n";

    int testCases;
    fin >> testCases;

    for (int tc = 1; tc <= testCases; tc++) {
        int n;
        fin >> n;

        vector<vector<double>> A(n, vector<double>(n));
        vector<double> B(n);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++)
                fin >> A[i][j];
            fin >> B[i];
        }

        fout << "Test Case " << tc << " (" << n << "x" << n << ")\n";
        fout << "Matrix:\n";

        vector<vector<double>> aug(n, vector<double>(n + 1));
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++)
                aug[i][j] = A[i][j];
            aug[i][n] = B[i];
        }
        printMatrix(fout, aug, n, n + 1);

        vector<vector<double>> L(n, vector<double>(n, 0.0));
        vector<vector<double>> U(n, vector<double>(n, 0.0));

        bool singular = false;

        for (int i = 0; i < n; i++) {
            L[i][i] = 1.0;

            for (int j = i; j < n; j++) {
                double sum = 0;
                for (int k = 0; k < i; k++)
                    sum += L[i][k] * U[k][j];
                U[i][j] = A[i][j] - sum;
            }

            if (fabs(U[i][i]) < eps) {
                singular = true;
                break;
            }

            for (int j = i + 1; j < n; j++) {
                double sum = 0;
                for (int k = 0; k < i; k++)
                    sum += L[j][k] * U[k][i];
                L[j][i] = (A[j][i] - sum) / U[i][i];
            }
        }

        fout << "L matrix:\n";
        printMatrix(fout, L, n, n);

        fout << "U matrix:\n";
        printMatrix(fout, U, n, n);

        if (singular) {
            int rankA = 0, rankAug = 0;

            for (int i = 0; i < n; i++) {
                bool nzA = false, nzAug = false;

                for (int j = 0; j < n; j++)
                    if (fabs(aug[i][j]) > eps) nzA = true;
                if (fabs(aug[i][n]) > eps) nzAug = true;

                if (nzA) rankA++;
                if (nzA || nzAug) rankAug++;
            }

            fout << "Matrix A is singular (det A = 0)\n";
            if (rankA < rankAug)
                fout << "-> System has \"NO\" solution.\n\n";
            else
                fout << "-> System has infinite solutions.\n\n";
            continue;
        }

        vector<double> y(n);
        for (int i = 0; i < n; i++) {
            y[i] = B[i];
            for (int j = 0; j < i; j++)
                y[i] -= L[i][j] * y[j];
        }

        vector<double> x(n);
        for (int i = n - 1; i >= 0; i--) {
            x[i] = y[i];
            for (int j = i + 1; j < n; j++)
                x[i] -= U[i][j] * x[j];
            x[i] /= U[i][i];
        }

        fout << "UNIQUE Solution!\n";
        for (int i = 0; i < n; i++) {
            fout << "x" << i + 1 << " = "
                 << fixed << setprecision(3)
                 << x[i] << endl;
        }
        fout << endl;
    }

    fin.close();
    fout.close();
    return 0;
}
