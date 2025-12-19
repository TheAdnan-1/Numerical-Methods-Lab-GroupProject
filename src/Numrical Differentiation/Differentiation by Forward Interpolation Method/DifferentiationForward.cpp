#include <bits/stdc++.h>
using namespace std;

double trueDerivative(vector<double> &coeffs, double x) {
    double result = 0;
    int n = coeffs.size();
    for (int i = 1; i < n; i++) {
        result += i * coeffs[i] * pow(x, i - 1);
    }
    return result;
}

int main() {
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    int testCases;
    fin >> testCases;

    for (int tc = 1; tc <= testCases; tc++) {
        int n, deg;
        fin >> deg >> n;
        vector<double> coeffs(deg + 1);
        for (int i = 0; i <= deg; i++)
            fin >> coeffs[i];

        vector<double> x(n), y(n);
        for (int i = 0; i < n; i++)
            fin >> x[i];
        for (int i = 0; i < n; i++) {
            double val = 0;
            for (int j = 0; j <= deg; j++)
                val += coeffs[j] * pow(x[i], j);
            y[i] = val;
        }

        double h = x[1] - x[0];
        double diffPoint;
        fin >> diffPoint;

        vector<vector<double>> F(n, vector<double>(n));
        for (int i = 0; i < n; i++)
            F[i][0] = y[i];

        for (int j = 1; j < n; j++) {
            for (int i = 0; i < n - j; i++)
                F[i][j] = F[i+1][j-1] - F[i][j-1];
        }

        double u = (diffPoint - x[0]) / h;
        double derivative = F[0][1] 
                          + ((2*u - 1)/2.0) * F[0][2]
                          + ((3*u*u - 6*u + 2)/6.0) * F[0][3]
                          + ((4*u*u*u - 18*u*u + 22*u - 6)/24.0) * F[0][4]; // up to 4th order
                          
        double trueDeriv = trueDerivative(coeffs, diffPoint);
        double error = fabs((trueDeriv - derivative) / trueDeriv) * 100;

        fout << "Test Case " << tc << "\n";
        fout << "Polynomial: f(x) = ";
        for(int i=0;i<=deg;i++){
            fout << fixed << setprecision(4) << coeffs[i];
            if(i>0) fout << "x^" << i;
            if(i<deg) fout << " + ";
        }
        fout << "\nNumber of points: " << n << "\n";
        fout << "x-values: ";
        for(double xi:x) fout << fixed << setprecision(6) << xi << " ";
        fout << "\n";
        fout << "y-values: ";
        for(double yi:y) fout << fixed << setprecision(6) << yi << " ";
        fout << "\n";
        fout << "Step size (h): " << fixed << setprecision(6) << h << "\n";
        fout << "Differentiation point: " << fixed << setprecision(6) << diffPoint << "\n";

        fout << "Forward Difference Table:\n";
        for(int i=0;i<n;i++){
            fout << "Row " << i << ": ";
            for(int j=0;j<n-i;j++)
                fout << fixed << setprecision(6) << F[i][j] << " ";
            fout << "\n";
        }

        fout << "Approximate derivative : " << fixed << setprecision(6) << derivative << "\n";
        fout << "True derivative: " << fixed << setprecision(6) << trueDeriv << "\n";
        fout << "Percentage error: " << fixed << setprecision(6) << error << " %\n\n";
    }

    fin.close();
    fout.close();
    return 0;
}
