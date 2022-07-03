//Autor: Guilherme Machado   ---Projeto de HPC 2---

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<omp.h>

#define L 3000  //linhas da matriz
#define C 3000 //colunas da matriz
#define MAX 10000
#define ERRO 1e-15
#define R 400

double contorno1(double **U)
{	int i, j;
	
	
	for(i=1;i<L;i++)
	{	for(j=1;j<C;j++)
		{	// as bordas das placas são tudo zero
		   U[L-1][j] = 0; 
			U[0][j] = 0; 
			U[i][C-1] = 0;
			U[i][0] = 0;
		}
	}
}	

//termo de fonte
double P(double x)
{	double o;
	o=1.0;
	return (o);

}


double relaxacao(double **G, double **D)
{	int i, j, x = L/4, y = C/4;
	double  h = 0.0001, k = 0.0001, aux, a;
	

	for(i=1;i<L;i++)
	{	for(j=1;j<C;j++)
		{	//fontes circulares
		  if((pow(i-x,2) + pow(j-y,2)<=R*R))              //750     750
				G[i][j] =exp(-i*j);
				
			if((pow(i-2250,2) + pow(j-2250,2)<=R*R))          // 2250   2250
				G[i][j] =exp(-i*j);
			
			if((pow(i-750,2) + pow(j-2250,2)<=R*R))         // 750   2250 
				G[i][j] =exp(-i*j);
				
				
			if((pow(i-2250,2) + pow(j-750,2)<=R*R))          //2250  750  
				G[i][j] =exp(-i*j);
						
			if((pow(i-1500,2) + pow(j-1500,2)<=R*R))        // 1500 1500        
				G[i][j] =exp(-i*j);
				
		}		
	}					
                             
                     #pragma omp parallel for       
                     for(i=1;i<L-1;i++)
                     {   
		          
                       for(j=1;j<C-1;j++)
                       {
		
	                          aux=G[i][j];    
                            G[i][j] = (((G[i+1][j] + G[i-1][j])*k*k + (G[i][j-1] + G[i][j+1])*h*h - k*h*k*h*P(0))/(2.0*(h*h + k*k))); //discretização       
                            aux = fabs(aux - G[i][j]);
			     
	        	}

               	}
		if(aux>a)
			a = aux;
	return a;

}


void salva_arquivo(double **S)
{	FILE *p;
	int i, j;
	
	p = fopen("poisson2.dat", "w+");
	
	for(i=0;i<L;i++)
	{	for(j=0;j<C;j++)
			fprintf(p, "%lf\t", S[i][j]);
			
		fprintf(p,"\n");
	}
	
	fclose(p);
}

int main()
{	double **M, **A, erro = 1;
	int i, j, cont;
	
	//alocando as matrizes
	M = (double**)malloc(L*sizeof(double*));
	for(i=0;i<L;i++)	
		M[i] = (double*)malloc(C*sizeof(double*));
	
	A = (double**)malloc(L*sizeof(double*));
	
	for(i=0;i<C;i++)	
		A[i] = (double*)malloc(C*sizeof(double*));
		
	for(i=0;i<L;i++)
		for(j=0;j<C;j++)
			M[i][j] = A[i][j] = 0;

	contorno1(M);

        #pragma omp parallel shared(M,A,cont,erro) private(i,j)
	{

          	for(cont=1;erro>ERRO && cont<MAX;cont++)
          	{     
									erro = relaxacao(M, A);

	        }
       
        }
	//salvando os valores no arquivo
	salva_arquivo(M);



}

