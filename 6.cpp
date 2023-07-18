#include <bits/stdc++.h>
#define eps 1e-8
using namespace std;

double f(double x, double t)
{
    return 4 * x / t + pow(t, 4) * exp(t);
}

double x(double t)
{
    return pow(t, 4) * (exp(t) - exp(1));
}

double *AB4(double t0, double tn, double h, double *w)
{
    double *error;
    int n = (tn - t0) / h;
    error = new double[n + 1];

    double wi[4];
    wi[0] = x(t0);
    for (int i = 1; i < 4; i++)
        wi[i] = wi[i - 1] + h * f(wi[i - 1], t0 + i * h);

    double t_curr;
    for (int i = 4; i < n + 1; i++)
    {
        t_curr = t0 + (i)*h;
        *w = wi[3] + h * (55 * f(wi[3], t_curr - h) - 59 * f(wi[2], t_curr - 2 * h) + 37 * f(wi[1], t_curr - 3 * h) - 9 * f(wi[0], t_curr - 4 * h)) / 24;
        *w = wi[3] + h * (9 * f(*w, t_curr) + 19 * f(wi[3], t_curr - h) - 5 * f(wi[2], t_curr - 2 * h) + f(wi[1], t_curr - 3 * h)) / 24;

        for (int j = 0; j < 3; j++)
        {
            wi[j] = wi[j + 1];
        }
        // cout << *w << endl;
        wi[3] = *w;
        error[i] = *w - x(t_curr);
    }
    return error;
}
double *RK4(double t0, double tn, double h, double *w)
{
    double *error;
    int n = (tn - t0) / h;
    error = new double[n + 1];

    *w = x(t0);
    double t_curr, k1, k2, k3, k4;
    // *error[0] = 0;
    for (int i = 0; i < n; i++)
    {
        t_curr = i * h + t0;
        error[i] = x(t_curr) - *w;

        k1 = h * f(*w, t_curr);
        k2 = h * f(*w + k1 / 2, t_curr + h / 2);
        k3 = h * f(*w + k2 / 2, t_curr + h / 2);
        k4 = h * f(*w + k3, t_curr + h);
        *w += (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    }

    error[n] = x(t_curr + h) - *w;
    return error;
}

int main(int argc, char const *argv[])
{

    double *error_RK4;
    double t0 = 1, tn = 2, h = 0.05;
    double *error_AB4;
    double w_RK4, w_AB4;
    error_RK4 = RK4(t0, tn, h, &w_RK4);
    error_AB4 = AB4(t0, tn, h, &w_AB4);
    cout << "Absolute value of x at t = " << tn << " is " << x(2) << "." << endl;
    cout << "Value of x by RK4 at t =  " << tn << " is " << w_RK4 << "." << endl;
    cout << "Value of x by AB4 at t = " << tn << " is " << w_AB4 << "." << endl;
    cout << "t\t"
         << "\nAbsolute Error By RK4\t"
         << "Absolute Error by AB4\n";

    for (int i = 0; i <= (tn - t0) / h; i++)
    {
        cout << t0 + i * h << "\t" << error_RK4[i] << "\t\t\t" << error_AB4[i] << endl;
    }

    delete (error_AB4);
    delete (error_RK4);
    return 0;
}