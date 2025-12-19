#include <bits/stdc++.h>
using namespace std;
const double eps = 1e-9;

void printMatrix(ofstream &fout, vector<vector<double>> &mat, int n, int m) {
    fout << "Matrix:\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            fout << setw(10) << fixed << setprecision(3) << mat[i][j] << " ";
        }
        fout << endl;
    }
    fout << endl;
}

int main() {
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    fout << "GAUSS JORDAN ELIMINATION\n\n";

    int testCases;
    fin >> testCases;
    fout << "Number of test cases: " << testCases << "\n\n";

    for (int tc = 1; tc <= testCases; tc++) {

        int n;
        fin >> n;

        vector<vector<double>> mat(n, vector<double>(n + 1));

        for (int i = 0; i < n; i++) {
            for (int j = 0; j <= n; j++) {
                fin >> mat[i][j];
            }
        }

        fout << "\nTest Case " << tc << " (" << n << "x" << n << ")\n";
        printMatrix(fout, mat, n, n + 1);

        int row = 0;
        for (int col = 0; col < n && row < n; col++) {

            int pivot = row;
            for (int i = row; i < n; i++) {
                if (fabs(mat[i][col]) > fabs(mat[pivot][col]))
                    pivot = i;
            }

            if (fabs(mat[pivot][col]) < eps)
                continue;

            swap(mat[row], mat[pivot]);

            double pivotValue = mat[row][col];
            for (int j = 0; j <= n; j++) {
                mat[row][j] /= pivotValue;
            }

            for (int i = 0; i < n; i++) {
                if (i == row) continue;
                double factor = mat[i][col];
                for (int j = 0; j <= n; j++) {
                    mat[i][j] -= factor * mat[row][j];
                }
            }

            row++;
        }

        int rankA = 0, rankAug = 0;

        for (int i = 0; i < n; i++) {
            bool nonZeroA = false;
            bool nonZeroAug = false;

            for (int j = 0; j < n; j++) {
                if (fabs(mat[i][j]) > eps)
                    nonZeroA = true;
            }
            if (fabs(mat[i][n]) > eps)
                nonZeroAug = true;

            if (nonZeroA)
                rankA++;
            if (nonZeroA || nonZeroAug)
                rankAug++;
        }

        if (rankA < rankAug) {
            fout << "No Solution!\n";
        }
        else if (rankA < n) {
            fout << "Infinite Solutions!\n";
        }
        else {
            fout << "Unique Solution!\n";
            for (int i = 0; i < n; i++) {
                fout << "x" << i + 1 << " = "
                     << fixed << setprecision(3)
                     << mat[i][n] << endl;
            }
        }
    }

    fin.close();
    fout.close();
    return 0;
}
