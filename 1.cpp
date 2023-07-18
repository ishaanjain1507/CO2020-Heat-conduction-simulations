#include <bits/stdc++.h>
using namespace std;

double f(double y, double t)
{
    return -(1 + t + t * t) - (2 * t + 1) * y - y * y;
}

double y(double t)
{
    return -t - 1 / (1 + exp(t));
}

double y1(double t)
{
    return -1 - t;
}
double Euler1(double t0, double tn, double h)
{
    double w;
    w = y(t0);

    int n = (tn - t0) / h;
    double t_curr;
    for (int i = 0; i < n; i++)
    {
        t_curr = i * h + t0;
        w += h * f(w, t_curr);
        // cout << y(t_curr + h) - w << endl;
    }
    // cout << w << endl;
    return w;
}
double Euler2(double t0, double tn, double h)
{
    double w;
    w = y1(t0);

    int n = (tn - t0) / h;
    for (int i = 1; i <= n; i++)
    {
        w += h * f(w, (i - 1) * h + t0);
    }
    // cout << w << endl;
    return w;
}
int main()
{
    double t0 = 0, tn = 3, h = 0.05;
    cout << "for part (a):" << endl;
    cout << "\tExact value = " << y(tn) << endl;
    cout << "\tValue calculated by Euler's method = " << Euler1(t0, tn, h) << " {when step size is 0.05.}" << endl;
    cout << "\tPercentage error = " << (y(tn) - Euler1(t0, tn, h)) * 100 / y(tn) << "%\n" << endl;

    cout << "for part (b):" << endl;
    cout << "\tExact value = " << y1(tn) << endl;
    cout << "\tValue calculated by Euler's method = " << Euler2(t0, tn, h) << " {when step size is 0.05.}" << endl;
    cout << "\tPercentage error = " << (y1(tn) - Euler2(t0, tn, h)) * 100 / y(tn) << "%" << endl;
    return 0;
}
