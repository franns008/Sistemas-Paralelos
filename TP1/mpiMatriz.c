#include <stdio.h>
#include <stdlib.h> 
#include <sys/time.h>
#include <stdbool.h>
#include <omp.h>

double maxGlobA,minGlobA,promGlobA,maxGlobB,minGlobB,promGlobB;
double resultado = 0.0;
int blocksize = 128;

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

void initMatrizTranspuestaMPI(double *matTrans, int n, double *mat){
    int i, j;
    #pragma omp for private(i,j)
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            matTrans[j*n + i] = mat[i*n + j];
        }
    }
}

void calcularMatrizPorNumeroMPI(double res, double *mat, int n){
    int i, j;
    #pragma omp for 
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            int index= i*n+j;
            mat[index] = mat[index] * res;
        }
    }
}

void multiplicarMatricesBloquesColumnasMPI(double *A, double *B, double *C, int N, int tamBloque) {
    int i, j, k, ii, jj, kk;
    #pragma omp for collapse(2) private(i,j,k,ii,jj,kk)
    for (i = 0; i < N; i += tamBloque) {
        for (j = 0; j < N; j += tamBloque) {
            for (k = 0; k < N; k += tamBloque) {
                for (ii = i; ii < i + tamBloque; ii++) {
                    for (jj = j; jj < j + tamBloque; jj++) {
                        double sum = 0.0;
                        int indexA = ii * N;
                        int indexB = jj * N;
                        for (kk = k; kk < k + tamBloque; kk++) {
                            sum += A[indexA + kk] * B[indexB+ kk];
                        }
                        C[indexA+ jj] += sum;
                    }
                }
            }
        }
    }
}


void calcularSumaMatricesMPI(double *a, double *b, double *res, int n){
    int i, j;
    #pragma omp for private(i,j)
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            int index= i*n+j;
            res[index] = a[index] + b[index];
        }
    }
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
    /*
    Si las matrices A, B, y C son definidas como matrices con todos sus elementos siendo 1,
    la ecuaciÃ³n deja una matriz de NxN con todos sus elementos siendo N
    */
   for (int i = 0; i < n * n; i++) {
        if (mat[i] != n) {
            printf("Error en la validacion de la matriz en posiciÃ³n %d\n", i);
            return;
        }
    }
    printf("El calculo de la ecuacion es correcto\n");
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

int main(int argc, char *argv[]) {
     // Tamaño del bloque para la multiplicación de matrices
    
    if (argc != 3){
        printf("Uso: %s <N> <numThreads>\n", argv[0]);
        return 1;
    }
    int N = atoi(argv[1]);
    int numThreads = atoi(argv[2]);

    if(N <= 0 || numThreads <= 0) {
        printf("Los valores de N y numThreads deben ser mayores a 0.\n");
        return 1;
    }

    double *A, *B, *C, *bTrans, *res1, *res2, *R;
    double maxGlobA, minGlobA, promGlobA, maxGlobB, minGlobB, promGlobB;
    double resultado = 0.0;
    double maxB = -1.0, minB = 9999999.0, sumB = 0.0;
    double maxA = -1.0, minA = 9999999.0, sumA = 0.0;

    A = (double*)malloc(sizeof(double) * N * N);
    B = (double*)malloc(sizeof(double) * N * N);
    C = (double*)malloc(sizeof(double) * N * N);
    bTrans = (double*)malloc(sizeof(double) * N * N);
    res1 = (double*)malloc(sizeof(double) * N * N);
    res2 = (double*)malloc(sizeof(double) * N * N);
    R = (double*)malloc(sizeof(double) * N * N);

    initValorMatrizFila(A, N, 1);
    initValorMatrizColumna(B, N, 1);
    initValorMatrizFila(C, N, 1); 
    initValorMatrizColumna(res1, N, 0);
    initValorMatrizColumna(res2, N, 0);

    long base = N / numThreads;
    long extra = N % numThreads;

    double timetick = dwalltime();
    #pragma omp parallel num_threads(numThreads)
    {
        
    initMatrizTranspuestaMPI(bTrans, N, B);

    #pragma omp for reduction(max:maxA) reduction(min:minA) reduction(+:sumA)
    for (int i = 0; i < N * N; i++) {
        double val = A[i];
        if (val > maxA) maxA = val;
        if (val < minA) minA = val;
        sumA += val;
    }

    #pragma omp for reduction(max:maxB) reduction(min:minB) reduction(+:sumB)
    for (int i = 0; i < N * N; i++) {
        double val = B[i];
        if (val > maxB) maxB = val;
        if (val < minB) minB = val;
        sumB += val;
    }
    
    double promGlobA = sumA / (N * N);
    double promGlobB = sumB / (N * N);
    
    multiplicarMatricesBloquesColumnasMPI(A, B, res1, N, blocksize);
    multiplicarMatricesBloquesColumnasMPI(C, bTrans, res2, N, blocksize);
    
    resultado = ((maxA * maxGlobB) - (minGlobA * minGlobB)) / (promGlobA * promGlobB);
    
    calcularMatrizPorNumeroMPI(resultado, res1, N);
    calcularSumaMatricesMPI(res1, res2, R, N);
    
    }
    timetick = dwalltime() - timetick;
    printf("Tiempo en segundos %f \n", timetick);
    printf("valor N y blocksize: %d %d\n", N, blocksize);
    printf("Cantidad de threads: %d\n", numThreads);
    validar(R, N);   

    

    return 0;
}