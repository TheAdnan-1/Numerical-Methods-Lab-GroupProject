#include <bits/stdc++.h>
using namespace std;

int degree;
double coeff[20];

// Polynomial function
double f(double x)
{
    double result = 0;
    for (int i = 0; i <= degree; i++)
    {
        result += coeff[i] * pow(x, degree - i);
    }
    return result;
}

int main()
{
    int t;
    cin >> t; // number of test cases

    while (t--)
    {
        double x0, x1, x2, tol;
        int maxIter, iter = 0;

        // Input polynomial
        cin >> degree;
        for (int i = 0; i <= degree; i++)
        {
            cin >> coeff[i];
        }

        // Input parameters
        cin >> x0 >> x1;
        cin >> tol;
        cin >> maxIter;

        bool found = false;

        while (iter < maxIter)
        {
            if (f(x1) - f(x0) == 0)
            {
                break;
            }

            // Secant formula
            x2 = x1 - (f(x1) * (x1 - x0)) / (f(x1) - f(x0));

            if (fabs(x2 - x1) < tol)
            {
                cout << "Root: " << x2 << endl;
                found = true;
                break;
            }

            x0 = x1;
            x1 = x2;
            iter++;
        }

        if (!found)
        {
            cout << "Not Converged" << endl;
        }
    }

    return 0;
}
