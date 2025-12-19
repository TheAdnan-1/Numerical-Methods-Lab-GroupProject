#include <bits/stdc++.h>
using namespace std;

double f(double x, double y) {
    return x*y + y;
}

int main() {
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    fout << "RUNGE KUTTA 4th ORDER METHOD\n\n";
    fout << "Function: f(x, y) = x*y + y\n\n";

    int T;
    fin >> T;
    fout << "Total Test Cases: " << T << "\n\n";

    for(int tc = 1; tc <= T; tc++) {
        double x0, y0, xn, h;
        fin >> x0 >> y0 >> xn >> h;

        int steps = (int)((xn - x0) / h);
        double x = x0;
        double y = y0;

        fout << "----------------------------------\n";
        fout << "TEST CASE #" << tc << "\n\n";
        fout << fixed << setprecision(3);
        fout << "Initial x0: " << x0 << ", y0: " << y0 << "\n";
        fout << "Final x: " << xn << "\n";
        fout << "Step size (h): " << h << "\n";
        fout << "Number of steps: " << steps << "\n\n";

        for(int i = 0; i < steps; i++) {
            double k1 = h * f(x, y);
            double k2 = h * f(x + h/2.0, y + k1/2.0);
            double k3 = h * f(x + h/2.0, y + k2/2.0);
            double k4 = h * f(x + h, y + k3);

            y = y + (k1 + 2*k2 + 2*k3 + k4)/6.0;
            x = x + h;
        }

        fout << "Result:\n";
        fout << "y(" << xn << ") = " << fixed << setprecision(3) << y << "\n\n";
    }

    fin.close();
    fout.close();
    return 0;
}
