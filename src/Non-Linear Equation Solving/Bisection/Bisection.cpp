#include <bits/stdc++.h>
using namespace std;

double f(double x)
{
    return x * x * x - x - 2;
}

int main()
{
    double a, b, mid;
    int n;

    cin >> a >> b; // input interval
    cin >> n;      // number of iterations

    if (f(a) * f(b) >= 0)
    {
        cout << "Invalid interval";
        return 0;
    }

    for (int i = 0; i < n; i++)
    {
        mid = (a + b) / 2;

        if (f(mid) == 0)
            break;
        else if (f(a) * f(mid) < 0)
            b = mid;
        else
            a = mid;
    }

    cout << mid;
    return 0;
}
