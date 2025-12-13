#include <bits/stdc++.h>
using namespace std;
#define arr vector<double>
double fun(const arr& coeff, double x) {
    double result = 0;
    int n = coeff.size(),i;
    for ( i = 0; i < n; i++) result += coeff[i] * pow(x, n - i - 1);
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
    for ( i = 0; i <= degree; i++) fin >> coeff[i];
    fin >> n;
    arr x(n), y(n);
    for ( i = 0; i < n; i++){
        fin >> x[i];
        y[i] = fun(coeff, x[i]);
    }
    double value;
    fin >> value;
    vector<arr> div(n, arr(n));
    for (i = 0; i < n; i++) div[i][0] = y[i];
    for (j = 1; j < n; j++) {
        for (i = 0; i < n - j; i++) {
            div[i][j] = (div[i + 1][j - 1] - div[i][j - 1]) /
                        (x[i + j] - x[i]);
        }
    }
    double result = div[0][0];
    double term = 1.0;
    for (i = 1; i < n; i++) {
        term *= (value - x[i - 1]);
        result += term * div[0][i];
    }
    fout << "Interpolated value at x = " << value << " is " << result << endl;
    }
    fin.close();
    fout.close();
    return 0;
}
