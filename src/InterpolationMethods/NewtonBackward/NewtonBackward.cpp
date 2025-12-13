#include <bits/stdc++.h>
using namespace std;
#define arr vector<double>
double fun(const arr& coeff, double x) {
    double result = 0;
    int n = coeff.size(),i;
    for (i = 0; i < n; i++) {
        result += coeff[i] * pow(x, n - i - 1);
    }
    return result;
}

int main() {
    ifstream fin("input.txt");
    ofstream fout("output.txt");
    int t;
    fin>>t;
    while(t--){
    int degree,n,i,j;
    fin >> degree;
    arr coeff(degree + 1);
    for (i = 0; i <= degree; i++) fin >> coeff[i];
    fin >> n;
    arr x(n), y(n);
    for ( i = 0; i < n; i++) {
        fin >> x[i];
        y[i] = fun(coeff, x[i]);
    }
    double value;
    fin >> value;
    vector<arr> diff(n, vector<double>(n));
    for ( i = 0; i < n; i++) diff[i][0] = y[i];
    for ( j = 1; j < n; j++) {
        for (int i = n - 1; i >= j; i--) {
            diff[i][j] = diff[i][j - 1] - diff[i - 1][j - 1];
        }
    }
    double h = x[1] - x[0];
    double u = (value - x[n - 1]) / h;
    double result = diff[n - 1][0];
    double u_term = 1;
    double fact = 1;
    for ( i = 1; i < n; i++) {
        u_term *= (u + (i - 1));
        fact *= i;
        result += (u_term * diff[n - 1][i]) / fact;
    }
    fout << "Interpolated value at x = " << value
         << " is " << result << endl;
    }
    fin.close();
    fout.close();
    return 0;
}
