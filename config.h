#pragma once
#include <math.h>
struct{
    unsigned N=5,M=5; //points count
    double x0 = 0;
    double xN = M_PI;
    double t0 = 0;
    double tM = M_PI;
    double tau = M_PI/2;
    double a = 1;
    double s = 0.8;
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
