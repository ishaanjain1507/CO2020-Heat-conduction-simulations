#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void exact(int nx, double *x, double *T_exact, double t_curr)
{
    for (int i = 0; i < nx; i++)
    {
        T_exact[i] = erf((x[i] - 0.5) / (2 * sqrt(t_curr)));
    }
}

void grid(int nx, double xst, double xen, double *x, double *dx)
{
    int i;

    *dx = (xen - xst) / (double)(nx - 1);

    for (i = 0; i < nx; i++)
        x[i] = (double)i * (*dx); // ensure x[0] == 0.0 and x[nx-1] == 1.0
}

void enforce_bcs(int nx, double *x, double *T)
{
    T[0] = -1.0;
    T[nx - 1] = 1.0;
}

void set_initial_condition(int nx, double *x, double *T)
{
    int i;

    for (i = 0; i < nx; i++)
        T[i] = tanh((x[i] - 0.5) / 0.05);

    enforce_bcs(nx, x, T); // ensure BCs are satisfied at t = 0
}

void get_rhs(int nx, double dx, double *x, double *T, double *rhs)
{
    int i;
    double dxsq = dx * dx;

    // compute rhs. For this problem, d2T/dx2
    rhs[1] = (0.91666667 * (T[0] - T[1]) + 0.5 * (T[2] - T[1]) + 0.333333 * (T[3] - T[1]) - 0.0833333 * (T[4] - T[1])) / dxsq;

    for (i = 2; i < nx - 2; i++)
    {
        rhs[i] = (-T[i + 2] + 16 * T[i + 1] - 30 * T[i] + 16 * T[i - 1] - T[i - 2]) / (12 * dxsq);
        // printf("%f ", rhs[i]);
    }
    // printf("\n");
    rhs[nx - 2] = (0.91666667 * (T[nx - 1] - T[nx - 2]) + 0.5 * (T[nx - 3] - T[nx - 2]) + 0.3333333 * (T[nx - 4] - T[nx - 2]) - 0.08333333 * (T[nx - 5] - T[nx - 2])) / dxsq;

    // for (int i = 0; i < nx; i++)
    // {
    //     printf("%f ", T[i]);
    // }
    // printf("\n");
}

void timestep_Euler(int nx, double dt, double dx, double *x, double *T, double *rhs)
{

    int i;

    // compute rhs
    get_rhs(nx, dx, x, T, rhs);

    // (Forward) Euler scheme
    for (i = 0; i < nx; i++)
        T[i] = T[i] + dt * rhs[i]; // T^(it+1)[i] = T^(it)[i] + dt * rhs[i];

    // set Dirichlet BCs
    enforce_bcs(nx, x, T);
}

void output_soln(int nx, int it, double tcurr, double *x, double *T, double *T_exact)
{
    int i;
    FILE *fp;
    char fname[100];

    sprintf(fname, "T_x_%04d.dat", it);
    // printf("\n%s\n", fname);

    fp = fopen(fname, "w");
    for (i = 0; i < nx; i++)
        fprintf(fp, "%lf %lf %lf\n", x[i], T[i], T_exact[i]);
    fclose(fp);

    // printf("Done writing solution for time step = %d\n", it);
}

double get_l2err_norm(int nx, double *T, double *T_exact)
{
    double l2err = 0.0, val;
    int i, j;

    for (i = 0; i < nx; i++)
    {
        val = T[i] - T_exact[i];
        l2err += val * val;
    }
    // printf("l2err = %lf\n", l2err);
    l2err = l2err / (double)(nx);
    l2err = sqrt(l2err);

    return l2err;
}

int main()
{

    int nx;
    double *x, *T, *rhs, tst, ten, xst, xen, dx, dt, tcurr, *T_exact;
    int i, it, num_time_steps, it_print;
    FILE *fp, *fp1;

    // read inputs
    fp = fopen("input.in", "r");
    fscanf(fp, "%d\n", &nx);
    fscanf(fp, "%lf %lf\n", &xst, &xen);
    fscanf(fp, "%lf %lf\n", &tst, &ten);
    fclose(fp);

    fp1 = fopen("./data.dat", "a+");
    printf("Inputs are: %d %lf %lf %lf %lf\n", nx, xst, xen, tst, ten);

    x = (double *)malloc(nx * sizeof(double));
    T = (double *)malloc(nx * sizeof(double));
    T_exact = (double *)malloc(nx * sizeof(double));

    rhs = (double *)malloc(nx * sizeof(double));
    for (int j = 10; j < nx; j = j+50)
    {
        grid(j, xst, xen, x, &dx); // initialize the grid
        // exact(j, xst, xen, x, T_exact);
        set_initial_condition(j, x, T); // initial condition

        // prepare for time loop
        dt = 3.125e-6;                                 // Ensure r satisfies the stability condition
        num_time_steps = (int)((ten - tst) / dt) + 1; // why add 1 to this?
        it_print = num_time_steps / 10;               // write out approximately 10 intermediate results

        // start time stepping loop
        for (it = 0; it < num_time_steps; it++)
        {
            tcurr = tst + (double)it * dt;

            exact(j, x, T_exact, tcurr);
            timestep_Euler(j, dt, dx, x, T, rhs); // update T

            // output soln every it_print time steps
            // if (it % it_print == 0)
            // output_soln(j, it, tcurr, x, T, T_exact);
        }
        double norm;
        norm = get_l2err_norm(j, T, T_exact);

        fprintf(fp1, "%d %lf\n", j, norm);
        printf("%d", j);
    }
    fclose(fp1);
    free(rhs);
    free(T);
    free(x);

    return 0;
}
