#include<stdio.h>
#include<stdlib.h> 
#include<sys/time.h>
#include <stdbool.h>

double *A,*B,*C,*bTrans,*res1,*res2,*R;

void initValorMatrizFila(double *mat, int n, double val)
{
  int i, j;      
	  for (i = 0; i < n; i++)
	  {
		for (j = 0; j < n; j++)
		{
		  mat[i*n + j] = val;
		}
	  }
}
void initValorMatrizColumna(double *mat, int n, double val)
{
  int i, j;      
	  for (i = 0; i < n; i++)
	  {
		for (j = 0; j < n; j++)
		{
		  mat[j*n + i] = val;
		}
	  }
}

void initValorMatrizTranspuesta(double *matTrans, int n, double *mat){
    int i, j;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            matTrans[j*n + i] = mat[i*n + j];
        }
    }
}

// Calcula el maximo, minimo y promedio de la matriz
void calcularMaximoMinimoPromedio(double *mat, int n, double *max, double *min, double *prom){
    int i, j;
    *max = -1;
    *min = 999999999;   
    *prom = 0;  
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            //se calcula el maximo, minimo y promedio de la matriz
            if (mat[i*n + j] > *max)
            {
                *max = mat[i*n + j];
            }
            if (mat[i*n + j] < *min)
            {
                *min = mat[i*n + j];
            }
            *prom += mat[i*n + j];
        }
    }
    *prom = *prom / (n*n);
}

void multiplicarMatricesBloquesColumnasB(double *A, double *B, double *C, int N, int tamBloque) {
    // Moverme entre fila de subbloques en la matriz A
    for (int i = 0; i < N; i += tamBloque) {

        // Me selecciona la proxima columna de "subbloques" de B, y A queda igual en la misma fila de "subbloques"
        for (int j = 0; j < N; j += tamBloque) {


            // Moverme al subbloque de la derecha en A, y en B al subloque de abajo en B
            for (int k = 0; k < N; k += tamBloque) {


                // selecciona la fila de la submatriz de la izquierda
                for (int ii = i; ii < i + tamBloque; ii++) {


                    // selecciona la columna de la submatriz derecha, manteniendo la misma fila de la submatriz izquierda
                    for (int jj = j; jj < j + tamBloque; jj++) {
                        double aux = 0.0; // Variable auxiliar para acumular la suma



                        // Este for va iterando entre las filas de la SUBmatriz izquierda y columnas de la SUBmatriz derecha
                        for (int kk = k; kk < k + tamBloque; kk++) {
                            aux += A[ii N + kk] * B[jj * N + kk];
                        }
                        C[ii * N + jj] += aux; // Asignar el resultado a C
                    }
                }
            }
        }
    }
}

void sumarMatrices(double *a, double *b, double *res, int n){
    int i, j;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            res[i*n + j] = a[i*n + j] + b[i*n + j];
        }
    }
}

void calcularMatrizPorNumero(double res, double *mat, int n){
    int i, j;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            mat[i*n + j] = mat[i*n + j] * res;
        }
    }
}

void printMatriz(double *mat, int n, bool transpuesta){
    int i, j;
    if(transpuesta){
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                printf("%f ", mat[j*n + i]);
            }
            printf("\n");
            
        }
    }else{
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                printf("%f ", mat[i*n + j]);
            }
            printf("\n");
            
        }
    }
    printf("\n");
    printf("\n");
}


double dwalltime(){
    double sec;
    struct timeval tv;

    gettimeofday(&tv,NULL);
    sec = tv.tv_sec + tv.tv_usec/1000000.0;
    return sec;
}

void validar(double *mat, int n){
    int i;
    
    for (i = 0; i < n; i++)
    {
        if (mat[i]!=n){
            printf("Error en la validacion de la matriz\n");
            printf("Error en la posicion %d\n", i);
            return;
        }
    }
    printf("El calculo de la ecuacion es correcto\n");
}


int main(int argc, char *argv[]) {
    
    if(argc != 2) {
        printf("Error: Se esperaba un argumento, el tamaño de la matriz y el blocksize \n");
    }

    int N = atoi(argv[1]);
    int blocksize = atoi(argv[2]);
    if(N <= 0 || blocksize <= 0) {
        printf("Error: El tamaño de la matriz y el blocksize deben ser mayores a 0\n");
        return -1;
    }
    
    A=(double*)malloc(sizeof(double)*N*N);
    B=(double*)malloc(sizeof(double)*N*N);
    C=(double*)malloc(sizeof(double)*N*N);
    
    bTrans=(double*)malloc(sizeof(double)*N*N);
    res1=(double*)malloc(sizeof(double)*N*N);
    res2=(double*)malloc(sizeof(double)*N*N);
    R=(double*)malloc(sizeof(double)*N*N);
    
    initValorMatrizFila(A, N, 1);
    initValorMatrizColumna(B, N, 1);
    initValorMatrizFila(C, N, 1); 
    initValorMatrizColumna(res1, N, 0);
    initValorMatrizColumna(res2, N, 0);

    double maxA, minA, promA, maxB, minB, promB;
    
    
    double timetick = dwalltime();
    
    initValorMatrizTranspuesta(bTrans, N,B);
    
    calcularMaximoMinimoPromedio(A, N, &maxA, &minA, &promA);
    calcularMaximoMinimoPromedio(B, N, &maxB, &minB, &promB);
    
    multiplicarMatricesBloquesColumnasB(A, B,res1, N,blocksize);
    multiplicarMatricesBloquesColumnasB(C, bTrans,res2, N, blocksize);
    

    double resultado = ((maxA*maxB)-(minA*minB))/(promA*promB);
    calcularMatrizPorNumero(resultado, res1, N);
    sumarMatrices(res2, res1, R, N);  

    printf("Tiempo en segundos %f\n", dwalltime() - timetick);
    printf("valor N y blocksize: %d %d\n", N, blocksize);

    validar(R, N);

    free(A);
    free(B);
    free(C);
    free(bTrans);
    free(res1);
    free(res2);
    free(R);
    
    return 0;
}