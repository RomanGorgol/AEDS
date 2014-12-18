#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

unsigned N=5,M=11; //points count
double x0 = 0,xN = M_PI, t0 = 0,tM = M_PI,  tau = M_PI/2 ,a = 1;
double h,d;
double s = 0.8; //вес между j и j+1 слоем. Если >0.5, то схема гарантированно устойчива
//k = floor(tau/d);%сколько раз шаг укладывается в величину запаздывания

double * g_U; //локальный результат

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

unsigned index(unsigned x,unsigned t) 
{
	return t*N + x;
//	return x*M+t;
}
double t(int i){
	return i*d+t0;
}
double x(int i){
	return i*h+x0;
}

int main(int argc, char* argv[]) {
	if (argc != 5) {
//		fprintf(stderr, "Usage: %s N input_file iterations output_file\n", argv[0]);
//		return 1;
	}
	int rank, size;
	//config

	// initialize globals
	h=(xN-x0)/(N-1);
	d=(tM-t0)/(M-1);
	
	unsigned k = floor(tau/d);

	g_U=(double*)malloc(sizeof(double)*N*M);
	
	unsigned i,j;
	for(i=0; i<M; ++i)
	{
		g_U[index(0,i)]=border_fun(t0+i*d);
	}

	for(i=0; i<N; ++i)
	{
		g_U[index(i,0)]=init_fun(x0+i*h,t0);
	}
	double lambda = tau/d - d;

	double * ut_j = (double*)malloc(sizeof(double)*N);
	double * ut_j1 = (double*)malloc(sizeof(double)*N);

	for(j=0; j<M-1; ++j)
	{
		if(j>=k+1)
		{
			ut_j1[0] = border_fun( t(j+1)-tau);
			ut_j[0]  = border_fun(t(j)-tau);
			for( i = 1; i<N; ++i)
			{
				ut_j1[i] = lambda*g_U[index(i,j+1-k)] + (1-lambda)*g_U[index(i,j-k)];
				ut_j[i]  = lambda*g_U[index(i,j-k)]   + (1-lambda)*g_U[index(i,j-1-k)];
			}
		}
		else if(j == k)
			{
				ut_j1[0] = border_fun(t(j+1)-tau);
				ut_j[0]  = init_fun( x(0), t(j)-tau);
				for( i = 1; i<N; ++i)
				{
					ut_j1[i] = lambda*g_U[index(i,j+1-k)] + (1-lambda)*g_U[index(i,j-k)]; //%lambda*U(i,2) + (1-lambda)*U(i,1);
					ut_j[i]  = init_fun(x(i), t(j)-tau);
				}
			}
			else
			{
				for(i = 0; i<N; ++i)
				{
					ut_j1[i] = init_fun(x(i), t(j+1)-tau);
					ut_j[i] =  init_fun(x(i), t(j)-tau);
				}
			} 
		
		double f_1_j1 = heter_fun(x0, t(j+1), g_U[index(0,j+1)], ut_j1[0]);
		double dot_g_j1 = der_border_fun(t(j+1));
    
		double f_1_j = heter_fun( x0, t(j), g_U[index(0,j)], ut_j[0]);
		double dot_g_j = der_border_fun(t(j));
    
		double F_2_j = 0.5 * (heter_fun(x(0), t(j+1), g_U[index(0,j+1)], ut_j1[0]) 
			+ heter_fun(x(1), t(j), g_U[index(1,j)], ut_j[1]));
    
		g_U[index(1,j+1)] = d*h/(h + 2*a*s*d)* (g_U[index(1,j)]/d -
			a*s		/ (2*h)*(-4*g_U[index(0,j+1)] - 2*h/a*(f_1_j1 - dot_g_j1)) - 
			a*(1-s) / (2*h)*(-4*g_U[index(0,j)]	  - 2*h/a*(f_1_j - dot_g_j) + 4*g_U[index(1,j)]) +  F_2_j);

		for(i=2; i<N; i++)
		{
			double F_i_j = 0.5 * (heter_fun(x(i-1), t(j+1), g_U[index(i-1,j+1)], ut_j1[i-1]) + 
				heter_fun(x(i), t(j), g_U[index(i,j)], ut_j[i]));

			g_U[index(i,j+1)] = 2*d*h/(2*h+3*a*s*d)* (g_U[index(i,j)]/d -
				a*s/    (2*h)*(g_U[index(i-2,j+1)] - 4*g_U[index(i-1,j+1)]) -
				a*(1-s)/(2*h)*(g_U[index(i-2,j)] - 4*g_U[index(i-1,j)] + 3*g_U[index(i,j)]) + F_i_j);
		}
	}
	
	FILE *output;
	output = fopen("output.txt","w");
	
	for(i=0; i<N; i++)
	{
		for(j=0; j<M; j++)
		{
			fprintf(output,"%lf\t",g_U[index(i,j)]);
		}
		fprintf(output,"\n");
	}
	fclose(output);
	printf("OK\n");
	return 0;
}
