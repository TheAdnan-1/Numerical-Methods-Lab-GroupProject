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

        // Simpson's 3/8 rule requires n multiple of 3
        if (n % 3 != 0)
        {
            cout << "Invalid n (must be multiple of 3)" << endl;
            continue;
        }

        h = (b - a) / n;
        sum = f(a) + f(b);

        for (int i = 1; i < n; i++)
        {
            double x = a + i * h;

            if (i % 3 == 0)
                sum += 2 * f(x);
            else
                sum += 3 * f(x);
        }

        double result = (3 * h / 8) * sum;
        cout << "Integral: " << result << endl;
    }

    return 0;
}
