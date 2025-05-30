#include <stdio.h>
#include <stdlib.h> 
#include <sys/time.h>
#include <stdbool.h>
#include <pthread.h>
#include <limits.h>

double *A,*B,*C,*bTrans,*res1,*res2,*R;
double maxGlobA, minGlobA, sumaGlobalA, maxGlobB, minGlobB, sumaGlobalB;

double resultado = 0.0;
int N, numThreads, tamBloque = 128;
long size, base, extra;

pthread_barrier_t barrier; // Barrera para sincronización de hilos
pthread_mutex_t mutexA, mutexB;


void initValorMatrizFila(double *mat, int n, double val){
  int i, j, indexI;      
	  for (i = 0; i < n; i++)
	  {
        indexI = i*n;
		for (j = 0; j < n; j++)
		{
		  mat[indexI + j] = val;
		}
	  }
}

void initValorMatrizColumna(double *mat, int n, double val){
  int i, j;      
	  for (i = 0; i < n; i++)
	  {
		for (j = 0; j < n; j++)
		{
		  mat[j*n + i] = val;
		}
	  }
}

void multiplicarMatricesBloquesColumnasBThread(double *A, double *B, double *C, int N, long inicio, long fin) {
    int indexII, indexJJ;    
    // Moverme entre fila de subbloques en la matriz A
    for (int i = inicio; i < fin; i += tamBloque) {

        // Me selecciona la proxima columna de "subbloques" de B, y A queda igual en la misma fila de "subbloques"
        for (int j = 0; j < N; j += tamBloque) {


            // Moverme al subbloque de la derecha en A, y en B al subloque de abajo en B
            for (int k = 0; k < N; k += tamBloque) {


                // selecciona la fila de la submatriz de la izquierda. Es importante el ii<fin para no salir del rango de la matriz 
                for (int ii = i; ii < i + tamBloque && ii < fin; ii++)  {
                    indexII = ii * N;

                    // selecciona la columna de la submatriz derecha, manteniendo la misma fila de la submatriz izquierda
                    for (int jj = j; jj < j + tamBloque; jj++) {
                        double aux = 0.0; // Variable auxiliar para acumular la suma
                        indexJJ = jj * N;
                        // Este for va iterando entre las filas de la SUBmatriz izquierda y columnas de la SUBmatriz derecha
                        for (int kk = k; kk < k + tamBloque; kk++) {
                            aux += A[indexII + kk] * B[indexJJ + kk];
                        }
                    
                        C[indexII + jj] += aux; // Asignar el resultado a C
                        
                    }
                }
            }
        }
    }
}

void sumarMatricesThread(double *a, double *b, double *res, int n, long inicio, long fin){
    int i, j, indexI, indexFinal;
    for (i = inicio; i < fin ; i++)
    {
        indexI = i*n;
        for (j = 0; j < n; j++)
        {
            indexFinal = indexI + j;
            res[indexFinal] = a[indexFinal] + b[indexFinal];
        }
    }
}

void initValorMatrizTranspuestaThread(double *matTrans, int n,double *mat,long inicio, long fin){

    int i, j, indexI;
    for (i = inicio; i < fin; i++)
    {
        indexI = i*n;
        for (j = 0; j < n; j++)
        {
            matTrans[j*n + i] = mat[indexI + j];
        }
    }
}

void calcularMaximoMinimoPromedioThread(double *mat, int n, double *max, double *min, double *prom,long inicio, long fin){
    int i, j, indexI;
    double elemento;
    *max = -1;
    *min = 999999999;   
    *prom = 0;  
    for (i = inicio; i < fin; i++)
    {
        indexI = i*n;
        for (j = 0; j < n; j++)
        {
            elemento = mat[indexI+j];
            //se calcula el maximo, minimo y promedio de la matriz
            if (elemento > *max)
            {
                *max = elemento;
            }
            if (elemento < *min)
            {
                *min = elemento;
            }
            
            *prom += elemento;
        }
    }
}

void calcularMatrizPorNumeroThread(double res, double *mat, int n, long inicio, long fin){
    int i, j, indexI, indexFinal;
    for (i = inicio; i < fin; i++)
    {
        indexI = i*n;
        for (j = 0; j < n; j++)
        {
            indexFinal = indexI + j;
            mat[indexFinal] = mat[indexFinal] * res;
        }
    }
}

void * pthreadResolucion(void *arg) {
    int thread_id = *(long *)arg;
    free(arg);
    long inicio = thread_id * base + (thread_id < extra ? thread_id : extra);
    long fin = inicio + base + (thread_id < extra ? 1 : 0);

    double maxA, minA, promA, maxB, minB, promB;
   

    initValorMatrizTranspuestaThread(bTrans, N,B,inicio,fin);

    calcularMaximoMinimoPromedioThread(A, N, &maxA, &minA, &promA,inicio,fin);
    calcularMaximoMinimoPromedioThread(B, N, &maxB, &minB, &promB,inicio,fin);
   
    pthread_mutex_lock(&mutexB); // Bloqueamos el mutex
    if (maxB > maxGlobB) {
        maxGlobB = maxB;
    }
    if (minB < minGlobB) {
        minGlobB = minB;
    }
    if (maxA > maxGlobA) {
        maxGlobA = maxA;
    }
    if (minA < minGlobA) {
        minGlobA = minA;
    }
    sumaGlobalA += promA;
    sumaGlobalB += promB;
    pthread_mutex_unlock(&mutexB); // Desbloqueamos el mutex

    multiplicarMatricesBloquesColumnasBThread(A, B,res1, N,inicio,fin);
    pthread_barrier_wait(&barrier); // Esperar a que todos los hilos terminen antes de continuar ==> Escalar y Bt
    multiplicarMatricesBloquesColumnasBThread(C, bTrans,res2, N ,inicio,fin);
    
    double promLocTotA = sumaGlobalA / (size);
    double promLocTotB = sumaGlobalB / (size);
    
    resultado = ((maxGlobA * maxGlobB) - (minGlobA * minGlobB)) / (promLocTotA * promLocTotB);

    // Suma y producto escalar
    
    calcularMatrizPorNumeroThread(resultado, res1, N,inicio,fin);
    sumarMatricesThread(res1, res2, R, N,inicio,fin);
    
    pthread_exit(NULL);
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
    /*
    Si las matrices A, B, y C son definidas como matrices con todos sus elementos siendo 1,
    la ecuación deja una matriz de NxN con todos sus elementos siendo N
    */
   for (int i = 0; i < n * n; i++) {
        if (mat[i] != n) {
            printf("Error en la validacion de la matriz en posicion %d\n", i);
            return;
        }
    }
    printf("El calculo de la ecuacion es correcto\n");
}


int main(int argc, char *argv[]) {
    if (argc != 3){
        printf("Uso: %s <N> <T>\n", argv[0]);
        return 1;
    }
    N = atoi(argv[1]);
    
    numThreads = atoi(argv[2]);
    if (N <= 0 || tamBloque <= 0) {
        printf("Los valores de N y tamBloque deben ser mayores que 0\n");
        return 1;
    }
    if (N % tamBloque != 0) {
        printf("El tamanio de la matriz debe ser divisible por el tamanio del bloque = %i \n", tamBloque);
        return 1;
    }
    if (N/tamBloque < numThreads) {
        tamBloque = 64;    
    }
    double maxA, minA, promA, maxB, minB, promB;
    size = N*N;
    pthread_t threads[numThreads]; 

    maxGlobA = maxGlobB = -99999999;
    minGlobA = minGlobB = 999999999;
    sumaGlobalA = sumaGlobalB = 0.0;

    A=(double*)malloc(sizeof(double)*size);
    B=(double*)malloc(sizeof(double)*size);
    C=(double*)malloc(sizeof(double)*size);
    
    bTrans=(double*)malloc(sizeof(double)*size);
    res1=(double*)malloc(sizeof(double)*size);
    res2=(double*)malloc(sizeof(double)*size);
    R=(double*)malloc(sizeof(double)*size);
    
    pthread_mutex_init(&mutexA, NULL); // Inicializamos el mutex
    pthread_mutex_init(&mutexB, NULL); // Inicializamos el mutex

    initValorMatrizFila(A, N, 1);
    initValorMatrizColumna(B, N, 1);
    initValorMatrizFila(C, N, 1); 
    initValorMatrizColumna(res1, N, 0);
    initValorMatrizColumna(res2, N, 0);
    pthread_barrier_init(&barrier, NULL, numThreads); // Inicializamos la barrera
    
    base = N / numThreads;    
    extra = N % numThreads;

    double timetick = dwalltime();

    for (int i = 0; i < numThreads; i++) {
        int *thread_id = malloc(sizeof(int));
        *thread_id = i;
        pthread_create(&threads[i], NULL, pthreadResolucion, (void *)thread_id);
    }

    
    for (int i = 0; i < numThreads; i++) {
        pthread_join(threads[i], NULL);
    }
    
    
    printf("Tiempo en segundos %f\n", dwalltime() - timetick);
    printf("valor N y blocksize: %d %d\n", N, tamBloque);
    printf("Threads: %d\n", numThreads);

    pthread_mutex_destroy(&mutexA); // Destruimos el mutex
    pthread_mutex_destroy(&mutexB); // Destruimos el mutex
    
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
