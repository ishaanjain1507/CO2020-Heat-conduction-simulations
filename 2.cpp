#include <bits/stdc++.h>
#define eps 10e-6
using namespace std;

double f(double x, double r, double k)
{
    return r * x * (1 - x / k) - x * x / (1 + x * x);
}

double euler(double w0, double t0, double h, double r, double k)
{
    double w = w0;
    double t_curr;
    int i = 0;
    while (1)
    {
        w += h * f(w, r, k);
        i++;
        if (abs(f(w, r, k)) < eps)
            break;
    }

    return w;
}

int main(int argc, char const *argv[])
{
    double w0 = 2.44, t0 = 0, h = 0.05, r = 0.4, k = 20;
    cout << "The population level became constant with the value " << euler(w0, t0, h, r, k) << " when initial population was " << w0 << endl;
    double w0_2 = 2.4;
    cout << "The population level became constant with the value " << euler(w0_2, t0, h, r, k) << " when initial population was " << w0_2 << endl;
    return 0;
}