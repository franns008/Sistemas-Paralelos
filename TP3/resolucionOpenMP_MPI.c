#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <sys/time.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>

#define MASTER 0

void calcularMatrizTranspuesta(double *B, double *Btrans, int n, int stripSize, int rank) {
    #pragma omp for schedule(static)
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < stripSize; j++) {
            Btrans[j * n + i] = B[i * n + (rank * stripSize + j)];
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
    MPI_Init(&argc, &argv);
    int i, j, n, rank, numProcs;
    double *A, *B, *C, *BtransLoc,*BtransTot, *res1, *res2, *R;
    double maxA, maxB, minA, minB, sumaA, sumaB;
    double localMaxA = -1.0, localMinA = 9999999.0, localSumaA = 0.0;
    double localMaxB = -1.0, localMinB = 9999999.0, localSumaB = 0.0;
    double promedioA, promedioB, escalar;
    double timetick, totalTime;
    int blockSize = 128;

    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(argc != 3){
        if(rank == MASTER) printf("Uso: %s <N> <numero de threads>\n", argv[0]);
        MPI_Finalize();
        return 1;
    }

    n = atoi(argv[1]);
    int numThreads = atoi(argv[2]);
    int stripSize = n / numProcs;

    if (rank == MASTER) {
        A = malloc(sizeof(double) * n * n);
        B = malloc(sizeof(double) * n * n);
        C = malloc(sizeof(double) * n * n);
        BtransLoc = malloc(sizeof(double) * n * n);
        BtransTot = malloc(sizeof(double) * n * n);
        res1 = calloc(n * n, sizeof(double));
        res2 = calloc(n * n, sizeof(double));
    } else {
        A = malloc(sizeof(double) * n * stripSize);
        B = malloc(sizeof(double) * n * n);
        C = malloc(sizeof(double) * n * stripSize);
        BtransLoc = malloc(sizeof(double) * n * stripSize);
        BtransTot = malloc(sizeof(double) * n * n);
        res1 = calloc(n * stripSize, sizeof(double));
        res2 = calloc(n * stripSize, sizeof(double));
    }
    R = malloc(sizeof(double) * n * n);

    if (rank == MASTER) {
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                A[i * n + j] = C[i * n + j] = B[i * n + j] = 1;
    }

    MPI_Scatter(A, n * stripSize, MPI_DOUBLE, A, n * stripSize, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Scatter(C, n * stripSize, MPI_DOUBLE, C, n * stripSize, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(B, n * n, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == MASTER) timetick = dwalltime();

    #pragma omp parallel num_threads(numThreads) shared(A, B, C, BtransTot, res1, res2, escalar)
    {
        calcularMatrizTranspuesta(B, BtransLoc, n, stripSize, rank);

        #pragma omp for reduction(max:localMaxA) reduction(min:localMinA) reduction(+:localSumaA)
        for (int i = 0; i < stripSize; i++) {
            int offsetI = i * n;
            for (int j = 0; j < n; j++) {
                double val = A[offsetI + j];
                if (val > localMaxA) localMaxA = val;
                if (val < localMinA) localMinA = val;
                localSumaA += val;
            }
        }

        #pragma omp for reduction(max:localMaxB) reduction(min:localMinB) reduction(+:localSumaB)
        for (int i = 0; i < stripSize; i++) {
            int offsetI = i * n;
            for (int j = 0; j < n; j++) {
                double val = B[offsetI + j];
                if (val > localMaxB) localMaxB = val;
                if (val < localMinB) localMinB = val;
                localSumaB += val;
            }
        }

        #pragma omp single
        {
            MPI_Gather(BtransLoc, n * stripSize, MPI_DOUBLE, BtransTot, n * stripSize, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
            MPI_Bcast(BtransTot, n * n, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
            MPI_Reduce(&localMaxA, &maxA, 1, MPI_DOUBLE, MPI_MAX, MASTER, MPI_COMM_WORLD);
            MPI_Reduce(&localMinA, &minA, 1, MPI_DOUBLE, MPI_MIN, MASTER, MPI_COMM_WORLD);
            MPI_Reduce(&localSumaA, &sumaA, 1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
            MPI_Reduce(&localMaxB, &maxB, 1, MPI_DOUBLE, MPI_MAX, MASTER, MPI_COMM_WORLD);
            MPI_Reduce(&localMinB, &minB, 1, MPI_DOUBLE, MPI_MIN, MASTER, MPI_COMM_WORLD);
            MPI_Reduce(&localSumaB, &sumaB, 1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);

            if (rank == MASTER) {
                promedioA = sumaA / (n * n);
                promedioB = sumaB / (n * n);
                escalar = ((maxA * maxB) - (minA * minB)) / (promedioA * promedioB);
            }

            MPI_Bcast(&escalar, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
        }

        multiplacionMatricesBloque(A, B, res1, blockSize, n, stripSize);
        multiplacionMatricesBloque(C, BtransTot, res2, blockSize, n, stripSize);
        multiplicarMatrizNumero(res1, escalar, stripSize, n);
        sumaMatrices(res1, res2, res1, n, stripSize);
    }

    MPI_Gather(res1, n * stripSize, MPI_DOUBLE, R, n * stripSize, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);

    if (rank == MASTER) {
        totalTime = dwalltime() - timetick;
        printf("El tiempo total es: %f\n", totalTime);

        for (int i = 0; i < n * n; i++) {
            if (R[i] != n) {
                printf("Error en el resultado, el valor es: %f en la posicion %i\n", R[i], i);
                printf("El resultado es incorrecto\n");
                break;
            }
        }
        printf("El resultado es correcto\n");
    }

    free(A); free(B); free(C); free(BtransLoc); free(BtransTot); free(res1); free(res2); free(R);

    MPI_Finalize();
    return 0;
}
