#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <sys/time.h>
#include <mpi.h>
#define MASTER 0

void calcularMaximoMinimoPromedio(double *A, int n, int stripSize, double *max, double *min, double *promedio){
    *promedio = 0;
    *max = -99999;
    *min = 99999;
    for (int i = 0; i < stripSize; i++) {
        int offsetI = i * n;
        for (int j = 0; j < n; j++) {
            double val = A[offsetI + j];
            if (val > *max) *max = val;
            if (val < *min) *min = val;
            *promedio += val;
        }
    }
}

// Cada proceso calcula su porción vertical de la transpuesta
void calcularMatrizTranspuesta(double *B, double *Btrans, int n, int stripSize, int rank) {
    int indexRank = rank * stripSize;
    int indexI;
    for (int i = 0; i < n; i++) {
        indexI = i * stripSize;
        for (int j = 0; j < stripSize; j++) {
            Btrans[indexI + j] = B[(indexRank + j) * n + i];
        }
    }
}



void sumaMatrices(double *A, double *B, double *C, int n, int stripSize){
    for (int i = 0; i < stripSize; i++) {
        int offsetI = i * n;
        for (int j = 0; j < n; j++) {
            C[offsetI + j] = A[offsetI + j] + B[offsetI + j];
        }
    }
}

void multiplacionMatricesBloque(double *A, double *B, double *C, int blockSize, int n, int stripSize){
    int offsetI, offsetJJ;
    double aux;
    for (int i = 0; i < stripSize; i += blockSize) {
        for (int j = 0; j < n; j += blockSize) {
            for (int k = 0; k < n; k += blockSize) {
                for (int ii = i; ii < i + blockSize && ii < stripSize; ii++) {
                    offsetI = ii * n;
                    for (int jj = j; jj < j + blockSize && jj < n; jj++) {
                        aux = 0.0;
                        offsetJJ = jj * n;
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
    MPI_Init(&argc, &argv);
    int i, j, n, rank, numProcs;
    double *A, *B, *C, *BtransLoc,*BtransTot, *res1, *res2, *R;
    double max[2], min[2], suma[2];
    double localMax[2] , localMin[2], localSuma[2];
    MPI_Status status;
    double commTime[8], commTimeMax[8], commTimeMin[8];

    int blockSize = 128;

    if ((argc != 2) || ((n = atoi(argv[1])) <= 0)) {
        printf("\nUsar: %s size \n  size: Dimension de la matriz y el vector\n", argv[0]);
        exit(1);
    }
    
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Request requests[9];
    int stripSize = n / numProcs;
    int size = n * n;
    int bufferStripSize = n * stripSize;
    if (stripSize < blockSize){
        blockSize = stripSize;
        if (rank == MASTER){
            printf("Se cambia el tamaño del bloque a %d para un mejor aprovechamiento de los procesos\n", blockSize);
        }
    }
    if (rank == MASTER) {
        printf("Numero de procesos: %d y tamaño de la matriz %i y blockSize %i \n", numProcs, n,blockSize);
    }
    if (rank == MASTER) {
        A = (double *)malloc(sizeof(double) * size);
        B = (double *)malloc(sizeof(double) * size);
        C = (double *)malloc(sizeof(double) * size);
        BtransLoc = (double *)malloc(sizeof(double) * size);
        BtransTot = (double *)malloc(sizeof(double) * size);
        res1 = (double *)calloc(size, sizeof(double));
        res2 = (double *)calloc(size, sizeof(double));
    } else {
        A = (double *)malloc(sizeof(double) * bufferStripSize);
        B = (double *)malloc(sizeof(double) * size);
        C = (double *)malloc(sizeof(double) * bufferStripSize);
        BtransLoc = (double *)malloc(sizeof(double) * bufferStripSize);
        BtransTot = (double *)malloc(sizeof(double) * size);
        res1 = (double *)calloc(bufferStripSize, sizeof(double));
        res2 = (double *)calloc(bufferStripSize, sizeof(double));
    }
    R = (double *)malloc(sizeof(double) * size);

    if (rank == MASTER) {
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                A[i * n + j] = 1;
                C[i * n + j] = 1;
            }
        }
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                B[i * n + j] = 1;
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    double tiempoEjecucion;
    if (rank == MASTER) {
        tiempoEjecucion = dwalltime();
    }

    commTime[0] = MPI_Wtime();
    MPI_Iscatter(A, bufferStripSize, MPI_DOUBLE, A, bufferStripSize, MPI_DOUBLE, MASTER, MPI_COMM_WORLD, &requests[0]);
    MPI_Ibcast(B, size, MPI_DOUBLE, MASTER, MPI_COMM_WORLD, &requests[1]);    
    MPI_Iscatter(C, bufferStripSize, MPI_DOUBLE, C, bufferStripSize, MPI_DOUBLE, MASTER, MPI_COMM_WORLD, &requests[2]); 
    
    commTime[1] = MPI_Wtime();
    
    MPI_Wait(&requests[0], MPI_STATUS_IGNORE);

    calcularMaximoMinimoPromedio(A, n, stripSize, &localMax[0], &localMin[0], &localSuma[0]);

    MPI_Wait(&requests[1], MPI_STATUS_IGNORE);
    calcularMaximoMinimoPromedio(B, n, stripSize, &localMax[1], &localMin[1], &localSuma[1]);
    calcularMatrizTranspuesta(B, BtransLoc, n, stripSize, rank);

    commTime[2] = MPI_Wtime();
    // Reducción de MaxA y MaxB    
    MPI_Ireduce(&localMax, &max, 2, MPI_DOUBLE, MPI_MAX, MASTER, MPI_COMM_WORLD, &requests[3]);
    
    // Reducción de MinA y MinB
    MPI_Ireduce(&localMin, &min, 2, MPI_DOUBLE, MPI_MIN, MASTER, MPI_COMM_WORLD, &requests[4]);


    // Reducción de sumaA y suma B
    MPI_Ireduce(&localSuma, &suma, 2, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD, &requests[5]);


    MPI_Iallgather(BtransLoc, bufferStripSize, MPI_DOUBLE, BtransTot, bufferStripSize, MPI_DOUBLE, MPI_COMM_WORLD, &requests[6]);
    commTime[3] = MPI_Wtime();

    
    multiplacionMatricesBloque(A, B, res1, blockSize, n, stripSize);
    
    MPI_Wait(&requests[2], MPI_STATUS_IGNORE);
    MPI_Wait(&requests[6], MPI_STATUS_IGNORE);
    multiplacionMatricesBloque(C, BtransTot, res2, blockSize, n, stripSize);
    
    double promedioA, promedioB, escalar;
    MPI_Wait(&requests[3], MPI_STATUS_IGNORE);
    MPI_Wait(&requests[4], MPI_STATUS_IGNORE);
    MPI_Wait(&requests[5], MPI_STATUS_IGNORE);
    
    if (rank == MASTER) {
    // ===================== FIN TERCERA ZONA =============================
        promedioA = suma[0] / size;
        promedioB = suma[1] / size;
        escalar = ((max[0] * max[1]) - (min[0] * min[1])) / (promedioA * promedioB);
    }

    commTime[4] = MPI_Wtime();
    MPI_Ibcast(&escalar, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD, &requests[7]);
    commTime[5] = MPI_Wtime();
    
    MPI_Wait(&requests[7], MPI_STATUS_IGNORE);
    
    multiplicarMatrizNumero(res1, escalar, stripSize, n);
    sumaMatrices(res1, res2, res1, n, stripSize);

    
    commTime[6] = MPI_Wtime();
    MPI_Igather(res1, bufferStripSize, MPI_DOUBLE, R, bufferStripSize, MPI_DOUBLE, MASTER, MPI_COMM_WORLD, &requests[8]);
    commTime[7] = MPI_Wtime();
    MPI_Wait(&requests[8], MPI_STATUS_IGNORE);
    if (rank == MASTER) {
        
        tiempoEjecucion = dwalltime()- tiempoEjecucion; 
    }
    double avgTime = (commTime[1] - commTime[0]) + (commTime[3] - commTime[2]) + (commTime[5] - commTime[4]) + (commTime[7] - commTime[6]);
    double totalTime = 0.0;
    printf("El tiempo de comunicacion del proceso %d es: %f\n", rank, avgTime);

    MPI_Reduce(&avgTime, &totalTime, 1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);

    if (rank == MASTER) {
        printf("El tiempo total es: %f\n", tiempoEjecucion);
        printf("El tiempo promedio de comunicacion es: %f\n", totalTime / numProcs);
        for (int i = 0; i < size; i++) {
            if (R[i] != n) {
                printf("Error en el resultado, el valor es: %f en la posicon %i\n", R[i], i);
                printf("El resultado es incorrecto\n");
                break;
            }
        }
        printf("El resultado es correcto\n");
    }

   if (rank == MASTER) {
        free(A);
        free(B);
        free(C);
        free(BtransLoc);
        free(BtransTot);
        free(res1);
        free(res2);
        free(R);
    } else {
        free(A);
        free(B);
        free(C);
        free(BtransLoc);
        free(BtransTot);
        free(res1);
        free(res2);
        free(R); 
    }

    MPI_Finalize();
    return 0;
}
