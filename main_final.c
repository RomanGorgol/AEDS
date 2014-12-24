#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <mpi.h>
#include "config.h"

/* MPI vars */

struct{
    int rank;
    int process_count;
    MPI_Request req[2];
} mpi;

int N,M; //points count
double x0,xN,t0,tM,tau,a,s = 0.8;
double h,d;
double * g_U; //локальный результат
double prev[2];    //граничные слева

//N - rows, M - columns
int index(int x,int t)
{
    return t*(N+2) + x + 2;
}
double t(int i){
    return i*d+t0;
}
double x(int i){
    return i*h+x0;
}


void configure(int argc, char* argv[]){
    /*
     * загрузка с хедера, возможна передача данных из аргументов или файла
     * писать с аргументов и файлов в config.*, т.к. занчения config.* используются также вне этой функции.
    */

    N=config.N;
    M=config.N;
    x0 = config.x0;
    xN = config.xN;
    t0 = config.t0;
    tM = config.tM;
    tau = config.tau;
    a = config.a;
    s = config.s;
    h=(xN-x0)/(N-1);
    d=(tM-t0)/(M-1);
}
void showMPIError(int err)
{
	char* errorText;
	int len;

	MPI_Error_string(err,errorText, &len);
}

int N_on_process(int process){
    int min_size=config.N / mpi.process_count;	//для удобства, минимальный размер разбиения у процесса

    int count_large=config.N % mpi.process_count;  //кол-во процессов, с разбиением, большем на 1.
    int first_large=mpi.process_count - count_large;

    return min_size+ ( (process < first_large )?0:1 );	//распределяем лоскуты поля, большие в _конец_, меньше нагрузка.
}

int start_from_on_process(int process)
{
    int min_size = config.N / mpi.process_count;

    int count_large = config.N % mpi.process_count;
    int first_large = mpi.process_count - count_large;

    int start = process*min_size;

    if( process > first_large){
        start+=process - first_large;
    }
    return start;
}

void configure_mpi(int argc, char* argv[]){
    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi.rank );
    MPI_Comm_size( MPI_COMM_WORLD, &mpi.process_count);
    N = N_on_process(mpi.rank);

    x0 = x(start_from_on_process(mpi.rank));
    xN = x0 + (N-1)*h;
}

/* начинаем прием данных, пишем в prev */
void start_recieving(){
    if(!mpi.rank) return;
	MPI_Status st;
    MPI_Irecv(prev, 2, MPI_DOUBLE, mpi.rank-1, MPI_ANY_TAG, MPI_COMM_WORLD,mpi.req);	//0ую получаем
}

void start_sending(double * what){
    if(mpi.rank<mpi.process_count-1){
		MPI_Status st;
        MPI_Isend(what, 2, MPI_DOUBLE, mpi.rank+1, MPI_ANY_TAG, MPI_COMM_WORLD,mpi.req+1);	//предпоследнюю посылаем
	}
}

void wait_read_request_completion()
{	
    MPI_Waitall(1, mpi.req, MPI_STATUSES_IGNORE);
}
/* ожидание завершения отправок */
void wait_request_completion()
{
    if(mpi.rank > 0 && mpi.rank < mpi.process_count-1)
	{
        MPI_Waitall(2, mpi.req, MPI_STATUSES_IGNORE);	
    }else if(mpi.rank==mpi.process_count-1)
	{
        wait_read_request_completion();
    }
	else
	{
        MPI_Waitall(1, mpi.req+1, MPI_STATUSES_IGNORE);
    }
}




int main(int argc, char* argv[]) {
	if (argc != 2) {
		fprintf(stderr, "Usage: %s output_file\n", argv[0]);
		return 1;
    }
	double startTime = MPI_Wtime();

	init_config();   
 	
	configure(argc,argv);    
	
    configure_mpi(argc,argv);


	
    int k = floor(tau/d);
    double lambda = tau/d - d;
    g_U=(double*)malloc(sizeof(double)*(N+2)*M);
    int i,j;
	
    if(mpi.rank > 0){
        //начинаем читать данные с предыдущей ноды
	    start_recieving();
    }
	else
	{
        for(i=0; i<M; ++i)
        {
            g_U[index(0,i)] = border_fun(t(i));
        }
    }
	
	
    i = 0;
    if(mpi.rank > 0)
	{
        i = -2;
    }


    for(; i<N; ++i)
    {
        g_U[index(i,0)] = init_fun(x(i),t0);
    }
	

    double * ut_j = (double*)malloc(sizeof(double)*(N+1));
    double * ut_j1 = (double*)malloc(sizeof(double)*(N+1));

    //чтобы можно пользовать [-1]
    ++ut_j;
    ++ut_j1;

    if(mpi.rank > 0)
	{
        wait_read_request_completion();
    }	
	
    for(j=0; j<M-1; ++j)
    {

        /* заполняем граничные и начинаем читать */
        if(mpi.rank > 0)
		{
            g_U[index(-2,j+1)] = prev[0];
            g_U[index(-1,j+1)] = prev[1];
        }
		
   		start_recieving();		
        
		if(j >= k+1)
        {
            /* получение края слева */
            i = -1;
            if(mpi.rank == 0)
            {
                i = 1;
                ut_j1[0] = border_fun( t(j+1)-tau);
                ut_j[0]  = border_fun(t(j)-tau);
            }
            for( ; i<N; ++i)
            {
                ut_j1[i] = lambda*g_U[index(i,j+1-k)] + (1-lambda)*g_U[index(i,j-k)];
                ut_j[i]  = lambda*g_U[index(i,j-k)]   + (1-lambda)*g_U[index(i,j-1-k)];
            }
        }
        else if(j == k)
            {
                i = -1;
                if(!mpi.rank)
                {
                    i = 1;
                    ut_j1[0] = border_fun(t(j+1)-tau);
                    ut_j[0]  = init_fun( x(0), t(j)-tau);
                }

                for(; i<N; ++i)
                {
                    ut_j1[i] = lambda*g_U[index(i,j+1-k)] + (1-lambda)*g_U[index(i,j-k)]; //%lambda*U(i,2) + (1-lambda)*U(i,1);
                    ut_j[i]  = init_fun(x(i), t(j)-tau);
                }
            }
            else
            {
                for(i =-1; i<N; ++i)
                {
                    ut_j1[i] = init_fun(x(i), t(j+1)-tau);
                    ut_j[i] =  init_fun(x(i), t(j)-tau);
                }
            }

        i=0;
        if(mpi.rank == 0)
		{
            i=2;
            double f_1_j1 = heter_fun(x0, t(j+1), g_U[index(0,j+1)], ut_j1[0]);
            double dot_g_j1 = der_border_fun(t(j+1));

            double f_1_j = heter_fun( x0, t(j), g_U[index(0,j)], ut_j[0]);
            double dot_g_j = der_border_fun(t(j));

            double F_2_j = 0.5 * (
                        heter_fun(x(0), t(j+1), g_U[index(0,j+1)],  ut_j1[0]) +
                        heter_fun(x(1), t(j),   g_U[index(1,j)],    ut_j[1])
                    );
			
            g_U[index(1,j+1)] = d*h/(h + 2*a*s*d)* (
                        g_U[index(1,j)]/d -
                        a*s		/ (2*h)*(-4*g_U[index(0,j+1)] - 2*h/a*(f_1_j1 - dot_g_j1)) -
                        a*(1-s) / (2*h)*(-4*g_U[index(0,j)]	  - 2*h/a*(f_1_j - dot_g_j) + 4*g_U[index(1,j)]) +
                        F_2_j
                    );			
        }
        for(; i<N; i++)
        {

			double F_i_j = 0.5 * (
                        heter_fun(x(i-1), t(j+1), g_U[index(i-1,j+1)], ut_j1[i-1]) +
                        heter_fun(x(i), t(j), g_U[index(i,j)], ut_j[i])
                    );

            g_U[index(i,j+1)] = 2*d*h/(2*h+3*a*s*d)* (g_U[index(i,j)]/d -
                a*s/    (2*h)*(g_U[index(i-2,j+1)] - 4*g_U[index(i-1,j+1)]) -
                a*(1-s)/(2*h)*(g_U[index(i-2,j)] - 4*g_U[index(i-1,j)] + 3*g_U[index(i,j)]) + F_i_j);
        }
		
		start_sending(&g_U[index(N-2,j+1)]);
		
		if (j < M-2)
		{
			wait_request_completion();
		}
    }
	double* U;
	int sizeU = (config.N+2*mpi.process_count)*config.M;
	if (mpi.rank == 0)
	{
		U = (double*)malloc(sizeof(double)*sizeU);
	}
	int* sizes = (int*)malloc(sizeof(int)*mpi.process_count);
	int* offsets = (int*)malloc(sizeof(int)*mpi.process_count);
	
	sizes[0] = M*(N_on_process(0)+2);
	offsets[0] = 0;
	for(i = 1; i < mpi.process_count; i++)
	{
		sizes[i] = M*(N_on_process(i)+2);
		offsets[i] = offsets[i-1]+sizes[i-1];		
	}

	
    MPI_Gatherv(g_U, M*(N+2), MPI_DOUBLE, U, sizes, offsets, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
	if (mpi.rank == 0)
	{
		printf("time: %lf\n",MPI_Wtime() - startTime);
    }  

	if (mpi.rank == 0)
	{
		FILE *output;
		output = fopen(argv[1],"w");
		
		int k = 0, offset = 0;
		for(k = 0; k < mpi.process_count; k++)
		{
			int width = N_on_process(k) + 2;
			for(i = 2; i<width; i++)
			{
				//(width+2)*M - количество эелментов, которое вернул к-ый процесс
				for(j = 0; j<M; ++j)
				{
					int ind = offset+j*width+i;
					fprintf(output,"%lf\t",U[ind]);
				}
				fprintf(output,"\n");
			}
			offset += width*M;
		}
		fclose(output);
	}

  //  printf("process %d, OK\n",mpi.rank);
  
	free(g_U);
    free(--ut_j);
    free(--ut_j1);
    
	MPI_Finalize();
	return 0;
}
