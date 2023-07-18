#include <bits/stdc++.h>
#define eps 1e-8
using namespace std;

double x(double t)
{
    return -1 / (t * log(2 * t));
}

double f(double x, double t)
{
    return x * x - x / t;
}

double *AB2(double t0, double tn, double h, double *w)
{
    double *error;
    int n = (tn - t0) / h;
    error = new double[n + 1];

    double wiplus1, wiminus1;
    wiminus1 = x(t0);
    *w = x(t0) + h * f(x(t0), t0);
    double t_curr, temp;
    error[0] = 0;
    for (int i = 0; i < n; i++)
    {
        t_curr = i * h + t0;
        wiplus1 = *w + h * (3 * f(*w, t_curr + h) - f(wiminus1, t_curr)) / 2;
        while (abs(wiplus1 - temp) < eps)
        {
            temp = wiplus1;
            wiplus1 = *w + h * (f(wiplus1, t_curr + 2 * h) + f(*w, t_curr + h)) / 2;
        }

        wiminus1 = *w;
        *w = wiplus1;
        error[i + 1] = *w - x(t_curr + h);
    }

    return error;
}
double *RK2(double t0, double tn, double h, double *w)
{
    double *error;
    int n = (tn - t0) / h;
    error = new double[n + 1];

    *w = x(t0);
    double t_curr, k1, k2;
    // *error[0] = 0;
    for (int i = 0; i < n; i++)
    {
        t_curr = i * h + t0;
        error[i] = x(t_curr) - *w;
        // cout << i << endl;

        k1 = h * f(*w, t_curr);
        k2 = h * f(*w + k1, t_curr + h);
        *w += (k1 + k2) / 2;
        // cout << i << endl;
    }

    error[n] = x(t_curr + h) - *w;
    return error;
}
int main(int argc, char const *argv[])
{
    double *error_RK2, *error_AB2;
    double t0 = 1, tn = 5, h = 0.05;
    double w_RK2, w_AB2;
    error_RK2 = RK2(t0, tn, h, &w_RK2);
    error_AB2 = AB2(t0, tn, h, &w_AB2);

    cout << "Absolute value of x at t = 5 is " << x(tn) << "." << endl;
    cout << "Value of x by RK2 at t = 5 is " << w_RK2 << "." << endl;
    cout << "Value of x by AB2 at t = 5 is " << w_AB2 << "." << endl;
    cout << "\nt\tAbsolute Error By RK2\tAbsolute Error by AB2\n";
    for (int i = 0; i <= (tn - t0) / h; i++)
    {
        cout << t0 + i*h << "\t" << error_AB2[i] << "\t\t\t" << error_RK2[i] << endl;
    }

    delete (error_AB2);
    delete (error_RK2);
    return 0;
}
