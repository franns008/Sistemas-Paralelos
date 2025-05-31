#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <sys/time.h>
#include <mpi.h>
#define MASTER 0

void calcularMaximoMinimoPromedio(double *A, int n, int stripSize, double *max, double *min, double *promedio){
    *promedio = 0;
    *max = A[0];
    *min = A[0];
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
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < stripSize; j++) {
            Btrans[i * stripSize + j] = B[(rank * stripSize + j) * n + i];
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
    double maxA, maxB, minA, minB, sumaA, sumaB;
    double localMaxA, localMinA, localSumaA;
    double localMaxB, localMinB, localSumaB;
    MPI_Status status;
    double timetick, totalTime;
    int blockSize = 128;

    if ((argc != 2) || ((n = atoi(argv[1])) <= 0)) {
        printf("\nUsar: %s size \n  size: Dimension de la matriz y el vector\n", argv[0]);
        exit(1);
    }
    
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    int stripSize = n / numProcs;

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
        A = (double *)malloc(sizeof(double) * n * n);
        B = (double *)malloc(sizeof(double) * n * n);
        C = (double *)malloc(sizeof(double) * n * n);
        BtransLoc = (double *)malloc(sizeof(double) * n * n);
        BtransTot = (double *)malloc(sizeof(double) * n * n);
        res1 = (double *)calloc(n * n, sizeof(double));
        res2 = (double *)calloc(n * n, sizeof(double));
    } else {
        A = (double *)malloc(sizeof(double) * n * stripSize);
        B = (double *)malloc(sizeof(double) * n * n);
        C = (double *)malloc(sizeof(double) * n * stripSize);
        BtransLoc = (double *)malloc(sizeof(double) * n * stripSize);
        BtransTot = (double *)malloc(sizeof(double) * n * n);
        res1 = (double *)calloc(n * stripSize, sizeof(double));
        res2 = (double *)calloc(n * stripSize, sizeof(double));
    }
    R = (double *)malloc(sizeof(double) * n * n);

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

    MPI_Scatter(A, n * stripSize, MPI_DOUBLE, A, n * stripSize, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Scatter(C, n * stripSize, MPI_DOUBLE, C, n * stripSize, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(B, n * n, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == MASTER) {
        timetick = dwalltime();
    }

    calcularMatrizTranspuesta(B, BtransLoc, n, stripSize, rank);
    MPI_Gather(BtransLoc, n * stripSize, MPI_DOUBLE, BtransTot, n * stripSize, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(BtransTot, n * n, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);

    calcularMaximoMinimoPromedio(A, n, stripSize, &localMaxA, &localMinA, &localSumaA);
    calcularMaximoMinimoPromedio(B, n, stripSize, &localMaxB, &localMinB, &localSumaB);
    multiplacionMatricesBloque(A, B, res1, blockSize, n, stripSize);
    multiplacionMatricesBloque(C, BtransTot, res2, blockSize, n, stripSize);

    MPI_Reduce(&localMaxA, &maxA, 1, MPI_DOUBLE, MPI_MAX, MASTER, MPI_COMM_WORLD);
    MPI_Reduce(&localMinA, &minA, 1, MPI_DOUBLE, MPI_MIN, MASTER, MPI_COMM_WORLD);
    MPI_Reduce(&localSumaA, &sumaA, 1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
    MPI_Reduce(&localMaxB, &maxB, 1, MPI_DOUBLE, MPI_MAX, MASTER, MPI_COMM_WORLD);
    MPI_Reduce(&localMinB, &minB, 1, MPI_DOUBLE, MPI_MIN, MASTER, MPI_COMM_WORLD);
    MPI_Reduce(&localSumaB, &sumaB, 1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);

    double promedioA, promedioB, escalar;

    if (rank == MASTER) {
        promedioA = sumaA / (n * n);
        promedioB = sumaB / (n * n);
        escalar = ((maxA * maxB) - (minA * minB)) / (promedioA * promedioB);
    }

    MPI_Bcast(&escalar, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    multiplicarMatrizNumero(res1, escalar, stripSize, n);
    sumaMatrices(res1, res2, res1, n, stripSize);

    MPI_Gather(res1, n * stripSize, MPI_DOUBLE, R, n * stripSize, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);

    if (rank == MASTER) {
        totalTime = dwalltime() - timetick;
        printf("El tiempo total es: %f\n", totalTime);

        for (int i = 0; i < n * n; i++) {
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
