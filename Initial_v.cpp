#include <bits/stdc++.h>
using namespace std;

double y(float t)
{
    return t * (1 + log(t));
}
double f(float y, float t)
{
    return 1 + y / t;
}
void Euler(double t0, double tn, double h)
{
    double wi_old, wi_new;
    wi_old = y(t0);

    int n = (tn - t0) / h;
    for (int i = 1; i <= n; i++)
    {
        wi_old += h * f(wi_old, (i - 1) * h + t0);
    }
    cout << wi_old << endl;
}

void RK2(double t0, double tn, double h)
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

    cout << w << endl;
}

void RK4(double t0, double tn, double h)
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

    cout << w << endl;
}
int main(int argc, char const *argv[])
{
    double t0 = 1, tn = 2, h = 0.05;
    cout << y(tn) << endl;
    Euler(t0, tn, h);
    RK2(t0, tn, h);
    RK4(t0, tn, h);
    return 0;
}
