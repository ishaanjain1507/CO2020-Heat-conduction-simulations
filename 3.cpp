#include <bits/stdc++.h>
#define sigma 5.67 * 10e-8
#define rho 7900
#define d 2 * 10e-3
#define eps 1e-6
using namespace std;

double CP(double T)
{
    return 0.162 * T + 446.47;
}

double f(double T, double Tf)
{
    return 2 * sigma * (pow(Tf, 4) - pow(T, 4)) / (rho * d * (CP(T) + 0.162 * T));
}

double RK4(double t0, double h, double T0, double Tf, double *T_eq)
{
    int i = 0;
    *T_eq = T0;
    double k1, k2, k3, k4;
    double t_curr;
    while (1)
    {
        t_curr = i * h + t0;
        k1 = h * f(*T_eq, Tf);
        k2 = h * f(*T_eq + k1 / 2, Tf);
        k3 = h * f(*T_eq + k2 / 2, Tf);
        k4 = h * f(*T_eq + k3, Tf);
        *T_eq += (k1 + 2 * k2 + 2 * k3 + k4) / 6;
        i++;

        if (abs(f(*T_eq, Tf)) < eps)
            break;
    }

    return i * h + t0;
}

int main(int argc, char const *argv[])
{
    double t0 = 0, h = 0.05, T0 = 300, Tf = 1500;
    double T_eq;
    cout << "In time " << RK4(t0, h, T0, Tf, &T_eq) << "seconds, the plate reaches equilibrium with the error of " << eps;
    cout << " in the slope of the T wrt time whem it reaches zero.\n"
         << "The final Temparature is ";
    cout << fixed << setprecision(8) << T_eq;
    cout << " K." << endl;
    return 0;
}
