#include <bits/stdc++.h>
using namespace std;
const double eps = 1e-9;

void getCofactor(vector<vector<double>> &A, vector<vector<double>> &temp, int p, int q, int n) {
    int i = 0, j = 0;
    for(int row = 0; row < n; row++) {
        for(int col = 0; col < n; col++) {
            if(row != p && col != q) {
                temp[i][j++] = A[row][col];
                if(j == n-1) { j = 0; i++; }
            }
        }
    }
}

double determinant(vector<vector<double>> &A, int n) {
    if(n == 1) return A[0][0];
    double det = 0;
    int sign = 1;
    vector<vector<double>> temp(n, vector<double>(n));
    for(int f = 0; f < n; f++) {
        getCofactor(A, temp, 0, f, n);
        det += sign * A[0][f] * determinant(temp, n-1);
        sign = -sign;
    }
    return det;
}

void adjoint(vector<vector<double>> &A, vector<vector<double>> &adj) {
    int n = A.size();
    if(n == 1) { adj[0][0] = 1; return; }
    int sign;
    vector<vector<double>> temp(n, vector<double>(n));
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            getCofactor(A, temp, i, j, n);
            sign = ((i+j)%2==0)? 1 : -1;
            adj[j][i] = sign * determinant(temp, n-1); 
        }
    }
}

vector<double> multiply(vector<vector<double>> &mat, vector<double> &v) {
    int n = mat.size();
    vector<double> res(n,0);
    for(int i=0;i<n;i++)
        for(int j=0;j<n;j++)
            res[i]+=mat[i][j]*v[j];
    return res;
}

void printMatrix(ofstream &fout, vector<vector<double>> &mat) {
    int n = mat.size(), m = mat[0].size();
    for(int i=0;i<n;i++){
        for(int j=0;j<m;j++)
            fout << setw(10) << fixed << setprecision(3) << mat[i][j] << " ";
        fout << endl;
    }
    fout << endl;
}

void printVector(ofstream &fout, vector<double> &v) {
    for(double x:v)
        fout << setw(10) << fixed << setprecision(3) << x << endl;
    fout << endl;
}

int main() {
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    fout << "MATRIX INVERSION METHOD\n\n";
    int testCases; 
    fin >> testCases;

    for(int tc=1;tc<=testCases;tc++){
        int n;
        fin >> n;
        vector<vector<double>> A(n, vector<double>(n));
        vector<double> B(n);

        fout << "Test Case " << tc << " (" << n << "x" << n << ")\n\n";

        fout << "Coefficient Matrix A:\n";
        for(int i=0;i<n;i++)
            for(int j=0;j<n;j++)
                fin >> A[i][j];
        printMatrix(fout, A);

        fout << "Vector b:\n";
        for(int i=0;i<n;i++) fin >> B[i];
        printVector(fout, B);

        double det = determinant(A,n);
        fout << "Determinant det(A): " << fixed << setprecision(3) << det << "\n\n";

        if(fabs(det)<eps){
            fout << "Matrix is singular (det=0)\n";
            fout << "-> System has no unique solution.\n\n";
            continue;
        }

        vector<vector<double>> adj(n, vector<double>(n));
        adjoint(A, adj);

        fout << "Adjoint Matrix adj(A):\n";
        printMatrix(fout, adj);

        vector<vector<double>> inv(n, vector<double>(n));
        for(int i=0;i<n;i++)
            for(int j=0;j<n;j++)
                inv[i][j] = adj[i][j]/det;

        fout << "Inverse Matrix A^-1:\n";
        printMatrix(fout, inv);

        vector<double> x = multiply(inv,B);
        fout << "Solution (x = A^-1 * b):\n";
        for(int i=0;i<n;i++)
            fout << "x" << i+1 << " = " << fixed << setprecision(3) << x[i] << endl;
        fout << endl;
    }

    fin.close();
    fout.close();
    return 0;
}
