#include <bits/stdc++.h>
using namespace std;

double f(double y, double t)
{
    return 7 * t * t - 4 * y / t;
}
double y(double t)
{
    return t * t * t + 1 / (t * t * t * t);
}
double RK2(double t0, double tn, double h)
{
    double k1, k2;
    double w;
    w = y(t0);
    int n = (tn - t0) / h;
    int t_curr;
    for (int i = 0; i < n; i++)
    {
        k1 = h * f(w, i * h + t0);
        k2 = h * f(w + k1, (i + 1) * h + t0);
        w += (k1 + k2) / 2;
    }

    return w;
}

double RK4(double t0, double tn, double h)
{
    double k1, k2, k3, k4;
    double w;
    w = y(t0);
    int n = (tn - t0) / h;
    int t_curr;
    for (int i = 0; i < n; i++)
    {
        k1 = h * f(w, i * h + t0);
        k2 = h * f(w + k1 / 2, i * h + t0 + h / 2);
        k3 = h * f(w + k2 / 2, i * h + t0 + h / 2);
        k4 = h * f(w + k3, i * h + t0 + h);
        w += (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    }

    return w;
}

int main()
{
    double t0 = 1, tn = 6, h = 0.05;

    cout << "Exact value = " << y(tn) << endl;
    cout << "\nRK2:" << endl;
    cout << "\tValue calculated by Euler's method = " << RK2(t0, tn, h) << " {when step size is 0.05.}" << endl;
    cout << "\tPercentage error = " << (y(tn) - RK2(t0, tn, h)) * 100 / y(tn) << "%" << endl;

    cout << "\nRK4:" << endl;
    cout << "\tValue calculated by Euler's method = " << RK4(t0, tn, h) << " {when step size is 0.05.}" << endl;
    cout << "\tPercentage error = " << (y(tn) - RK4(t0, tn, h)) * 100 / y(tn) << "%" << endl;
    return 0;
}
