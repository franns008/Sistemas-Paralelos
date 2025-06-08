#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <sys/time.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>

#define MASTER 0

void calcularMatrizTranspuesta(double *B, double *Btrans, int n, int stripSize, int rank) {
    int indexRank = rank * stripSize;
    int indexI;
    #pragma omp for schedule(static) private(indexI, indexRank)
    for (int i = 0; i < n; i++) {
        indexI = i * n;
        for (int j = 0; j < stripSize; j++) {
            Btrans[j * n + i] = B[indexI + (indexRank + j)];
        }
    }
}

void sumaMatrices(double *A, double *B, double *C, int n, int stripSize){
    #pragma omp for schedule(static)
    for (int i = 0; i < stripSize; i++) {
        int offsetI = i * n;
        for (int j = 0; j < n; j++) {
            C[offsetI + j] = A[offsetI + j] + B[offsetI + j];
        }
    }
}

void multiplacionMatricesBloque(double *A, double *B, double *C, int blockSize, int n, int stripSize){
    #pragma omp for schedule(static)
    for (int i = 0; i < stripSize; i += blockSize) {
        for (int j = 0; j < n; j += blockSize) {
            for (int k = 0; k < n; k += blockSize) {
                for (int ii = i; ii < i + blockSize && ii < stripSize; ii++) {
                    int offsetI = ii * n;
                    for (int jj = j; jj < j + blockSize && jj < n; jj++) {
                        double aux = 0.0;
                        int offsetJJ = jj * n;
                        for (int kk = k; kk < k + blockSize && kk < n; kk++) {
                            aux += A[offsetI + kk] * B[offsetJJ + kk];
                        }
                        C[offsetI + jj] += aux;
                    }
                }
            }
        }
    }
}

void multiplicarMatrizNumero(double *A, double escalar, int stripSize, int n){
    #pragma omp for schedule(static)
    for (int i = 0; i < stripSize; i++) {
        int indexI = i * n;
        for (int j = 0; j < n; j++) {
            A[indexI + j] *= escalar;
        }
    }
}

double dwalltime(){
    double sec;
    struct timeval tv;
    gettimeofday(&tv, NULL);
    sec = tv.tv_sec + tv.tv_usec / 1000000.0;
    return sec;
}

int main(int argc, char *argv[]){

    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    int i, j, n, rank, numProcs;
    double *A, *B, *C, *BtransLoc,*BtransTot, *res1, *res2, *R;
    double max[2], min[2], suma[2];
    double localMax[2] , localMin[2], localSuma[2];
    double commTimeMax[8];
    double commTimeMin[8];
    double commTime[8];
    double escalar;
    int blockSize = 128;
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    // ====================================VALIDACIONES===================================
    
    if(argc != 3){
        printf("Uso: %s <N> <numero de threads>\n", argv[0]);
        MPI_Finalize();
        return 1;
    }
    
    n = atoi(argv[1]);
    
    if (n <= 0){
        if (rank == MASTER) {
            printf("El valor de N debe ser mayor a 0. Se ingresó =>%i<=\n", n);
        }
        MPI_Finalize();
        return 1;    
    }
    int numThreads = atoi(argv[2]);
    if (numThreads <= 0){
        if (rank == MASTER) {
            printf("El valor de numThreads debe ser mayor a 0. Se ingresó =>%i<=\n", numThreads);
        }
        MPI_Finalize();
        return 1;    
    }
    
    int stripSize = n / numProcs;
   
    if ((stripSize / numThreads) < blockSize){
        if (rank == MASTER){
            printf("Cambiaremos el blocksize al valor %i para una ejecución balanceada \n", (stripSize / numThreads));
        }
        blockSize = (stripSize/ numThreads);
    }
    
    if (rank == MASTER) {
        printf("El tamaño de la matriz es: %i x %i\n", n, n);
        printf("El número de procesos por hilo es: %i\n", numProcs);
        printf("El número de threads es: %i\n", numThreads);
        printf("El tamaño del bloque es: %i\n", blockSize);
    }
    
    // ====================================FIN-VALIDACIONES===================================
    
    int size = n * n;
    int bufferSizeStrip = n * stripSize;
    
    if (rank == MASTER) {
        A = malloc(sizeof(double) * size);
        B = malloc(sizeof(double) * size);
        C = malloc(sizeof(double) * size);
        BtransLoc = malloc(sizeof(double) * size);
        BtransTot = malloc(sizeof(double) * size);
        res1 = calloc(size, sizeof(double));
        res2 = calloc(size, sizeof(double));
    } else {
        A = malloc(sizeof(double) * bufferSizeStrip);
        B = malloc(sizeof(double) * size);
        C = malloc(sizeof(double) * bufferSizeStrip);
        BtransLoc = malloc(sizeof(double) * bufferSizeStrip);
        BtransTot = malloc(sizeof(double) * size);
        res1 = calloc(bufferSizeStrip, sizeof(double));
        res2 = calloc(bufferSizeStrip, sizeof(double));
    }
    R = malloc(sizeof(double) * size);

    if (rank == MASTER) {
        for (i = 0; i < n; i++){
            int indexI = i * n;
            for (j = 0; j < n; j++){
                int indexFinal = indexI + j;
                A[indexFinal] = B[indexFinal] = C[indexFinal] = 1;
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    double tiempoEjecucion;
    if (rank == MASTER) {
        tiempoEjecucion = dwalltime();
    }
    commTime[0] = MPI_Wtime();

    // Esto no cuenta como parte del tiempo? <===========================================================================================================================
    MPI_Scatter(A, bufferSizeStrip, MPI_DOUBLE, A, bufferSizeStrip, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Scatter(C, bufferSizeStrip, MPI_DOUBLE, C, bufferSizeStrip, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(B, size, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    commTime[1] = MPI_Wtime();
    

    #pragma omp parallel num_threads(numThreads) shared(A, B, C, BtransTot, res1, res2, escalar)
    {
        calcularMatrizTranspuesta(B, BtransLoc, n, stripSize, rank);
        localMax[0] = localMax[1] = -1;
        localMin[0] = localMin[1] = 9999999;
        localSuma[0] = localSuma[1] = 0.0;

        #pragma omp for reduction(max:localMax[0]) reduction(min:localMin[0]) reduction(+:localSuma[0])
        for (int i = 0; i < stripSize; i++) {
            int offsetI = i * n;
            for (int j = 0; j < n; j++) {
                double val = A[offsetI + j];
                if (val > localMax[0]) localMax[0] = val;
                if (val < localMin[0]) localMin[0] = val;
                localSuma[0] += val;
            }
        }

        #pragma omp for reduction(max:localMax[1]) reduction(min:localMin[1]) reduction(+:localSuma[1])
        for (int i = 0; i < stripSize; i++) {
            int offsetI = i * n;
            for (int j = 0; j < n; j++) {
                double val = B[offsetI + j];
                if (val > localMax[1]) localMax[1] = val;
                if (val < localMin[1]) localMin[1] = val;
                localSuma[1] += val;
            }
        }

        #pragma omp single
        {
            commTime[2] = MPI_Wtime();
            MPI_Allgather(BtransLoc, n * stripSize, MPI_DOUBLE, BtransTot, n * stripSize, MPI_DOUBLE, MPI_COMM_WORLD);
            // Reducción a un valor del máximo, mínimo y sumas locales de A y B
            MPI_Reduce(&localMax, &max, 2, MPI_DOUBLE, MPI_MAX, MASTER, MPI_COMM_WORLD);
            MPI_Reduce(&localMin, &min, 2, MPI_DOUBLE, MPI_MIN, MASTER, MPI_COMM_WORLD);
            MPI_Reduce(&localSuma, &suma, 2, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
            commTime[3] = MPI_Wtime();


            // Con los datos exactos, realizar el calculo del escalar
            double promedioA, promedioB;

            if (rank == MASTER) {
                promedioA = suma[0] / (n * n);
                promedioB = suma[1] / (n * n);
                escalar = ((max[0] * max[1]) - (min[0] * min[1])) / (promedioA * promedioB);
            }
            
            // Informar a todos los procesos del escalar
            commTime[4] = MPI_Wtime();
            MPI_Bcast(&escalar, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
            commTime[5] = MPI_Wtime();
        }

        multiplacionMatricesBloque(A, B, res1, blockSize, n, stripSize);
        multiplacionMatricesBloque(C, BtransTot, res2, blockSize, n, stripSize);
        multiplicarMatrizNumero(res1, escalar, stripSize, n);
        sumaMatrices(res1, res2, res1, n, stripSize);
    }
    

    
    // Obtener el resultado final a partir de lo parcial de todos los nodos
    commTime[6] = MPI_Wtime();
    MPI_Gather(res1, bufferSizeStrip, MPI_DOUBLE, R, bufferSizeStrip, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    commTime[7] = MPI_Wtime();
    if(rank == MASTER) {
        tiempoEjecucion = dwalltime() - tiempoEjecucion;
    }    

    double avgTime = (commTime[1] - commTime[0]) + (commTime[7] - commTime[6]) +
                     (commTime[3] - commTime[2]) + (commTime[5] - commTime[4]);
    printf("Tiempo de comunicación del proceso %d: %f\n", rank, avgTime);
    double avgTimeTot = 0.0;
    MPI_Reduce(&avgTime, &avgTimeTot, 1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
    
    if (rank == MASTER) {
        avgTimeTot /= numProcs;
        printf("El tiempo total es: %f\n", tiempoEjecucion);
        printf("El tiempo promedio de comunicacion es: %f\n", avgTimeTot);
        for (int i = 0; i < n * n; i++) {
            if (R[i] != n) {
                printf("Error en el resultado, el valor es: %f en la posicion %i\n", R[i], i);
                printf("El resultado es incorrecto\n");
                break;
            }
        }
        printf("El resultado es correcto, con N= %i y numero de procesos = %i \n", n, numProcs);
    }

    free(A); free(B); free(C); free(BtransLoc); free(BtransTot); free(res1); free(res2); free(R);

    MPI_Finalize();
    return 0;
}
