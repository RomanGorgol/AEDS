#ifndef XONFIG_H
#define XONFIG_H
#include <math.h>

struct{
    int N;
	int M; //points count
    double x0;
    double xN;
    double t0;
    double tM;
    double tau;
    double a;
    double s;
} config;

double init_fun(double x,double t)
{
    return sin(x)*cos(t);
}
double border_fun(double t)
{
    return 0;
}
double der_border_fun(double t)
{
    return 0;
}
double heter_fun(double x,double t,double u,double ut_j)
{
    return cos(x)*cos(t)-ut_j;
}
void init_config()
{
    config.N = 1000000;
	config.M = 6000; //points count

    config.x0 = 0;
    config.xN = M_PI;
    config.t0 = 0;
    config.tM = M_PI;
    config.tau = M_PI/2;
    config.a = 1;
    config.s = 0.8;
}

#endif
