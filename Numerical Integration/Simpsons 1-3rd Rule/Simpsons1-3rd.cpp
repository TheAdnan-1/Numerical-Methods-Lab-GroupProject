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
    cin >> t;

    while (t--)
    {
        double a, b, h, sum;
        int n;

        // Input polynomial
        cin >> degree;
        for (int i = 0; i <= degree; i++)
        {
            cin >> coeff[i];
        }

        // Input limits and subintervals
        cin >> a >> b;
        cin >> n;

        // Simpson's 1/3 rule requires even n
        if (n % 2 != 0)
        {
            cout << "Invalid n (must be even)" << endl;
            continue;
        }

        h = (b - a) / n;
        sum = f(a) + f(b);

        for (int i = 1; i < n; i++)
        {
            double x = a + i * h;
            if (i % 2 == 0)
                sum += 2 * f(x);
            else
                sum += 4 * f(x);
        }

        double result = (h / 3) * sum;

        cout << "Integral: " << result << endl;
    }

    return 0;
}
