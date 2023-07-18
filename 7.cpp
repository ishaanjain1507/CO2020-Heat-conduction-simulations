#include <bits/stdc++.h>
using namespace std;

vector<double> f(vector<double> u)
{
    vector<double> a;
    a.push_back(u.at(1));
    a.push_back(-u.at(0) + (1 - u.at(0) * u.at(0)) * u.at(1));

    return a;
}
int main()
{
    vector<double> uo, u;
    uo = {0.5, 0.1};
    u = {0.5, 0.1};
    vector<double> w0, w1;
    double to, tn, h;
    int n;
    to = 0;
    tn = 30;
    h = 0.1;

    n = (int)((tn - to) / h);
    vector<double> temp;
    temp = {0, 0};
    w0.push_back(uo.at(0));
    w1.push_back(uo.at(1));
    vector<double> k1, k2;
    k1 = {0, 0};
    k2 = {0, 0};
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < 2; j++)
            k1[j] = h * f(u)[j];
        for (int j = 0; j < 2; j++)
            u[j] = u[j] + k1[j];
        for (int j = 0; j < 2; j++)
        {
            k2[j] = h * f(u)[j];
        }
        w0.push_back(w0.at(i) + (k1[0] + k2[0]) / 2);
        w1.push_back(w1.at(i) + (k1[1] + k2[1]) / 2);
    }

    cout << "t\tx" << endl;
    for (int i = 0; i < n + 1; i++)
    {
        cout << to + i * h << "\t" << w0[i] << endl;
    }

    cout << "The graph is plotted using the given data and the period is found using that graph by iterpolation. Pls refer report." << endl;
    return 0;
}