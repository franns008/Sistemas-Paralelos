#include<stdio.h>
#include<stdlib.h> 
#include<sys/time.h>
#include <stdbool.h>
#include <pthread.h>

double *A,*B,*C,*bTrans,*res1,*res2,*R;
double maxGlobA, minGlobA, promGlobA, maxGlobB, minGlobB, promGlobB;

// Mutex para sincronización
pthread_mutex_t mutexMaxA = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mutexMinA = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mutexPromA = PTHREAD_MUTEX_INITIALIZER;

pthread_mutex_t mutexMaxB = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mutexMinB = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mutexPromB = PTHREAD_MUTEX_INITIALIZER;

double escalar = 0.0;
int N, blocksize;
int numThreads;

typedef struct {
    double *matriz;
    int inicio, fin;
} datosHiloEscalar;

typedef struct {
  double *a, *b, *c;
  int start, end;
} ThreadArgs;

typedef struct {
    double *matriz;
    int inicio, fin;
    double escalar;
} datosEscalarPorMatriz;

typedef struct {
    int start;
    int end;
} argSumaMatrices;

typedef struct{
    int inicio, fin;
    double *original, *transpuesta;
} argTranspuesta;
/* Multiply square matrices, blocked version */
void *matmulblks(void *arg);

/* Multiply (block)submatrices */
void blkmul(double *ablk, double *bblk, double *cblk);

void initValorMatrizFila(double *mat, int n, double val){
  int i, j;      
	  for (i = 0; i < n; i++)
	  {
		for (j = 0; j < n; j++)
		{
		  mat[i*n + j] = val;
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


void sumarMatricesThread(double *a, double *b, double *res, int n, long inicio, long fin){
    int i, j;
    for (i = inicio; i < fin; i++)
    {
        for (j = 0; j < n; j++)
        {
            res[i*n + j] = a[i*n + j] + b[i*n + j];
        }
    }
}

void initValorMatrizTranspuestaThread(double *matTrans, int n,double *mat,long inicio, long fin){

    int i, j;
    for (i = inicio; i < fin; i++)
    {
        for (j = 0; j < n; j++)
        {
            matTrans[j*n + i] = mat[i*n + j];
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
    /*
    Si las matrices A, B, y C son definidas como matrices con todos sus elementos siendo 1,
    la ecuación deja una matriz de NxN con todos sus elementos siendo N
    */
   for ( i = 0; i < n * n; i++) {
        if (mat[i] != n) {
            printf("Error en la validacion de la matriz en posición %d\n", i);
            return;
        }
    }
    printf("El calculo de la ecuacion es correcto (dio todo %i)\n", n);
}


void *calculoEscalaresA(void *arg) {
    datosHiloEscalar *datos = (datosHiloEscalar *)arg;
    double *matriz = datos->matriz;
    int inicio = datos->inicio;
    int fin = datos->fin;
    double x;
    double maxLocal = -999999;
    double minLocal = 999999;
    double sumaParcial = 0;

     for (int i = inicio; i < fin; i++) {
        
        x = matriz[i];
        
        if (x > maxLocal){
            maxLocal = x;
        }
        if (x < minLocal){
            minLocal = x;
        }
        sumaParcial += x;
    }
    //-----------------------------------------------------
    // mutex de maxA
    pthread_mutex_lock(&mutexMaxA);
    if (maxLocal > maxGlobA){
        maxGlobA = maxLocal;
    }
    // demutex de maxA
    pthread_mutex_unlock(&mutexMaxA);
    //-----------------------------------------------------

    //-----------------------------------------------------
    // mutex de minA
    pthread_mutex_lock(&mutexMinA);
    if (minLocal < minGlobA) {
        minGlobA = minLocal;
    }
    // demutex de minA
    pthread_mutex_unlock(&mutexMinA);
    //-----------------------------------------------------

    //-----------------------------------------------------
    // mutex de promA
    pthread_mutex_lock(&mutexPromA);
    promGlobA += sumaParcial;
    // demutex de promA
    pthread_mutex_unlock(&mutexPromA);
    //-----------------------------------------------------
    pthread_exit(NULL);
}

void *calculoEscalaresB(void *arg) {
    datosHiloEscalar *datos = (datosHiloEscalar *)arg;
    double *matriz = datos->matriz;
    int inicio = datos->inicio;
    int fin = datos->fin;
    double x;
    double maxLocal = -999999;
    double minLocal = 999999;
    double sumaParcial = 0;

    for (int i = inicio; i < fin; i++) {
        
        x = matriz[i];
        
        if (x > maxLocal){
            maxLocal = x;
        }
        if (x < minLocal){
            minLocal = x;
        }
        sumaParcial += x;
    }
        
    
    
    //-----------------------------------------------------
    // mutex de maxB
    pthread_mutex_lock(&mutexMaxB);
    if (maxLocal > maxGlobB){
        maxGlobB = maxLocal;
    }
    // demutex de maxB
    pthread_mutex_unlock(&mutexMaxB);
    //-----------------------------------------------------

    //-----------------------------------------------------
    // mutex de minB
    pthread_mutex_lock(&mutexMinB);
    if (minLocal < minGlobB) {
        minGlobB = minLocal;
    }
    // demutex de minA
    pthread_mutex_unlock(&mutexMinB);
    //-----------------------------------------------------

    //-----------------------------------------------------
    // mutex de promB
    pthread_mutex_lock(&mutexPromB);
    promGlobB += sumaParcial;
    // demutex de promB
    pthread_mutex_unlock(&mutexPromB);
    //-----------------------------------------------------
    pthread_exit(NULL);
}


/* Multiply square matrices, blocked version */
void *matmulblks(void *arg)
{
    ThreadArgs *args = (ThreadArgs *)arg;
    int i, j, k;
    // Moverme entre fila de subbloques en la matriz A
    for (i = args->start; i < args->end; i += blocksize) {


        // Me selecciona la próxima columna de "subbloques" B, y A queda igual en la misma fila de "subbloques"
        for (j = 0; j < N; j += blocksize) {


            //Moverme al subbloque de la derecha en A, y en B al subbloque de abajo
            for (k = 0; k < N; k += blocksize) {
                blkmul(&args->a[i * N + k],
                       &args->b[j * N + k],
                       &args->c[i * N + j]);
            }
        }
    }

    pthread_exit(NULL);
}

void blkmul(double *ablk, double *bblk, double *cblk)
{
    int ii, jj, kk, indexII, indexJJ;
    for (ii = 0; ii < blocksize; ii++)
    {
        // Selecciona la columna de la submatriz derecha, dejando fija la fila de la submatriz izquierda
        for (jj = 0; jj < blocksize; jj++)
        {
              double cont = 0.0;
              indexII = ii*N;
              indexJJ = jj*N;
              // Este for va iterando por los elementos de la fila y la columna
              for  (kk = 0; kk < blocksize; kk++){
                cont += ablk[indexII + kk] * bblk[indexJJ + kk];
              } 
              cblk[indexII + jj] += cont;

        }
    }
}


void *transponerParalelo(void *arg) {
    argTranspuesta *args = (argTranspuesta *)arg;
    int inicio = args->inicio;
    int fin = args->fin;
    double *original = args->original;
    double *transpuesta = args->transpuesta;
    for (int i = inicio; i < fin; i ++){
        for (int j=0; j<N; j++){
            
            transpuesta[j*N + i] = original[i*N + j];        
        }    
    }
    pthread_exit(NULL);
}

void *escalarPorMatriz(void *arg) {
    datosEscalarPorMatriz *args = (datosEscalarPorMatriz *)arg;
    double *matriz = args->matriz;
    int inicio = args->inicio;
    int fin = args->fin;
    double escalar = args->escalar;

    for (int i = inicio; i < fin; i++) {
        for (int j = 0; j < N; j++) {
            matriz[i * N + j] *= escalar;
        }
    }


    pthread_exit(NULL);
}

void * sumaMatrizChunk(void * ptr){
    argSumaMatrices *args = (argSumaMatrices *) ptr;
    // args->start es como acceder a args.start, lo mismo que args->end seria args.end
    
    for (int i = args->start; i < args->end; i++) {
        for (int j=0; j < N; j ++){
            R[i*N + j] = res1[i*N + j] + res2[i*N + j];    
        }
    }
    pthread_exit(NULL);
}

int main(int argc, char *argv[]) {
    if (argc != 4){
        printf("Uso: %s <N> <blocksize> <T>\n", argv[0]);
        return 1;
    }
    N = atoi(argv[1]);
    blocksize = atoi(argv[2]);
    numThreads = atoi(argv[3]);
    if (N <= 0 || blocksize <= 0 || numThreads <= 0) {
        printf("Los valores de N, blocksize y numThreads deben ser mayores que 0\n");
        return 1;
    }
    if (N % blocksize != 0) {
        printf("El tamaño de la matriz debe ser divisible por el tamaño del bloque\n");
        return 2;
    }

    int size = N*N;
    pthread_t threads[numThreads]; 

    maxGlobA = maxGlobB = -99999999;
    minGlobA = minGlobB = 999999999;
    promGlobA = promGlobB = 0.0;
    int tamParte = size / numThreads;
    int filas = N / numThreads;
    
    A=(double*)malloc(sizeof(double)*size);
    B=(double*)malloc(sizeof(double)*size);
    C=(double*)malloc(sizeof(double)*size);
    
    bTrans=(double*)malloc(sizeof(double)*size);
    res1 = (double*)calloc(size, sizeof(double));

    res2 = (double*)calloc(size, sizeof(double));

    R = (double*)calloc(size, sizeof(double));

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Analizar si se podrian hacer en stack ==> numThreads siempre es chico < 16
    
    datosHiloEscalar *datos = (datosHiloEscalar *)malloc(sizeof(datosHiloEscalar) * numThreads);
    datosEscalarPorMatriz *datosEporM = (datosEscalarPorMatriz *)malloc(sizeof(datosEscalarPorMatriz) * numThreads);
    ThreadArgs *args = (ThreadArgs *)malloc(sizeof(ThreadArgs) * numThreads);
    argSumaMatrices *argSuma = (argSumaMatrices *)malloc(sizeof(argSumaMatrices) * numThreads);
    argTranspuesta *argTrans = (argTranspuesta *)malloc(sizeof(argTranspuesta) * numThreads);
    // Analizar si se podrian hacer en stack ==> numThreads siempre es chico < 16
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    initValorMatrizFila(A, N, 1);
    initValorMatrizColumna(B, N, 1);
    initValorMatrizFila(C, N, 1); 
    initValorMatrizColumna(res1, N, 0);
    initValorMatrizColumna(res2, N, 0);
    
    double timetick = dwalltime();
    
// .........................................................................................................
    // Calcular escalares en A
    for (int i = 0; i < numThreads; i++) {
        datos[i].matriz = A;
        datos[i].inicio = i * tamParte;

        // Si soy el ultimo hilo, agarrar el resto
        if (i == numThreads - 1) {
            datos[i].fin = size;        
        }
        else {
            datos[i].fin = (i + 1)* tamParte;
        }
        pthread_create(&threads[i], NULL, calculoEscalaresA, &datos[i]);
    }
    
    // Esperar a que los hilos de A terminen
    for (int i = 0; i < numThreads; i++) {
        pthread_join(threads[i], NULL);
    }


// .........................................................................................................
    // Calcular escalares en B
      for (int i = 0; i < numThreads; i++) {
            datos[i].matriz = B;
            datos[i].inicio = i * tamParte;

            // Si soy el ultimo hilo, agarrar el resto
            if (i == numThreads - 1) {
                datos[i].fin = size;        
            }
            else {
                datos[i].fin = (i + 1)* tamParte;
            }
            pthread_create(&threads[i], NULL, calculoEscalaresB, &datos[i]);
        }
        
        // Esperar a que los hilos de A terminen
    for (int i = 0; i < numThreads; i++) {
        pthread_join(threads[i], NULL);
    }
    
// ......................................................................................................... 
    promGlobA = promGlobA / size;
    promGlobB = promGlobB / size;
    if (promGlobA == 0 || promGlobB == 0) {
        fprintf(stderr, "Error: promedio cero, no se puede dividir.\n");
        exit(EXIT_FAILURE);
    }

    escalar = ((maxGlobA * maxGlobB) - (minGlobA * minGlobB)) / (promGlobA * promGlobB);
// ......................................................................................................... 

    // Multiplicar [A x B]
    int rows_per_thread = N / numThreads;
    
    for (int t = 0; t < numThreads; t++) {
        args[t].a = A;
        args[t].b = B;
        args[t].c = res1;
        args[t].start = t * rows_per_thread;
        
        if (t == numThreads-1){
            args[t].end = N;        
        }
        else {
            args[t].end = args[t].start + rows_per_thread;         
        }
        pthread_create(&threads[t], NULL, matmulblks, &args[t]);
    }
    
    for (int t = 0; t < numThreads; t++) {
        pthread_join(threads[t], NULL);
    }
    


    
// .........................................................................................................
    // Transponer B en paralelo
    // ESTO SE HACE POR FILAS; NO POR ELEMENTOS. Por eso NO SE DEBE USAR size ACA
    
    for (int t=0; t < numThreads; t++){
        
        argTrans[t].original = B;
        argTrans[t].transpuesta = bTrans;
        argTrans[t].inicio = t * filas;
        if (t == numThreads-1){
            argTrans[t].fin = N;        
        }
        else {
            argTrans[t].fin = argTrans[t].inicio + filas;         
        }
        pthread_create(&threads[t], NULL, transponerParalelo, &argTrans[t]);
    }

    for (int t = 0; t < numThreads; t++) {
        pthread_join(threads[t], NULL);
    }
    



// .........................................................................................................


    // Multiplicar [C x Bt]
    
    for (int t=0; t < numThreads; t++){
        args[t].a = C;
        args[t].b = bTrans;
        args[t].c = res2;
        args[t].start = t * rows_per_thread;
        
        if (t == numThreads-1){
            args[t].end = N;        
        }
        else {
            args[t].end = args[t].start + rows_per_thread;         
        }
        pthread_create(&threads[t], NULL, matmulblks, &args[t]);
    }
    
    for (int t = 0; t < numThreads; t++) {
        pthread_join(threads[t], NULL);
    }
// .........................................................................................................


    // Escalar x [res1]
    // SE DEBE DE HACER POR FILAS
    
    
    for (int t=0; t < numThreads; t++){
            datosEporM[t].matriz = res1;
            datosEporM[t].inicio = t * filas;
            // Si soy el ultimo hilo, agarrar el resto
            if (t == numThreads - 1) {
                datosEporM[t].fin = N;        
            }
            else {
                datosEporM[t].fin = datosEporM[t].inicio + filas;
            }
            datosEporM[t].escalar = escalar;
            pthread_create(&threads[t], NULL, escalarPorMatriz, &datosEporM[t]); 
    }
    
    for (int t = 0; t < numThreads; t++) {
        pthread_join(threads[t], NULL);
    }
    printf("...\n");
    
//.........................................................................................................

    // Sumar [res1] + [res2]
    // ESTO SE DEBE DE HACER POR FILAS
    
    for (int t = 0; t < numThreads; t++) {
            
        argSuma[t].start = t * filas;
        
        // si estoy en la ultima iteración
        if (t == numThreads - 1) {
            // El último thread toma hasta el final
            argSuma[t].end = N; 
        } else {
            // Si no agarro la parte que me corresponde
            argSuma[t].end = (t + 1) * filas;
        }

        pthread_create(&threads[t], NULL, sumaMatrizChunk, &argSuma[t]);
    }
    for (int i=0; i<numThreads; i++){
	    pthread_join(threads[i], NULL);
    }
    //.........................................................................................................
    
    printf("Tiempo en segundos %f\n", dwalltime() - timetick);
    printf("valor N y blocksize: %d %d\n", N, blocksize);
    printf("Numero de threads = %i \n", numThreads);
    validar(R,N);

    free(A); free(B); free(C); free(bTrans);
    free(res1); free(res2); free(R);
    free(datos); free(datosEporM); free(args);
    free(argSuma); free(argTrans);

    
    pthread_mutex_destroy(&mutexMaxA);
    pthread_mutex_destroy(&mutexMinA);
    pthread_mutex_destroy(&mutexPromA);
    pthread_mutex_destroy(&mutexMaxB);
    pthread_mutex_destroy(&mutexMinB);
    pthread_mutex_destroy(&mutexPromB);
    return 0;
}
