#include <bits/stdc++.h>
using namespace std;

// Evaluate polynomial f(x) given coefficients
double f(double x, double c[], int deg)
{
    double res = 0;
    for (int i = 0; i <= deg; i++)
        res += c[i] * pow(x, deg - i);
    return res;
}

int main()
{
    int T;
    cin >> T;

    cout << "Total Test Cases: " << T << "\n\n";

    for (int t = 1; t <= T; t++)
    {
        int deg;
        cin >> deg;

        double coef[20];
        for (int i = 0; i <= deg; i++)
            cin >> coef[i];

        double a = -10, b = 10;
        double tol = 0.001;
        int n = 50;
        double fa = f(a, coef, deg);
        double fb = f(b, coef, deg);

        // Print polynomial in readable form
        cout << "Function: f(x) = ";
        for (int i = 0; i <= deg; i++)
        {
            if (i > 0 && coef[i] >= 0)
                cout << "+";
            cout << coef[i];
            if (deg - i > 0)
                cout << "x^" << (deg - i) << " ";
        }
        cout << "\n";

        cout << "Degree: " << deg << "\n";
        cout << fixed << setprecision(6);
        cout << "Error Tolerance: " << tol << "\n";
        cout << "Search Interval: [" << a << ", " << b << "]\n";
        cout << "Root Approximation:\n";

        if (fa * fb > 0)
        {
            cout << "  Invalid interval (no sign change)\n";
        }
        else
        {
            double c;
            for (int i = 0; i < n; i++)
            {
                // False position formula
                c = (a * fb - b * fa) / (fb - fa);
                double fc = f(c, coef, deg);

                cout << "  Iter " << i + 1 << ": x = " << c << "\n";

                if (fabs(fc) < tol)
                    break;

                if (fa * fc < 0)
                {
                    b = c;
                    fb = fc;
                }
                else
                {
                    a = c;
                    fa = fc;
                }
            }
            cout << "  Final root ≈ " << c << "\n";
        }
    }
    return 0;
}