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
    int indexRank = rank * stripSize;
    int indexI;
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
    // ====================================VALIDACIONES===================================
    if(argc != 3){
        if(rank == MASTER) printf("Uso: %s <N> <numero de threads>\n", argv[0]);
        MPI_Finalize();
        return 1;
    }

    n = atoi(argv[1]);
    if (n <= 0){
        printf("El valor de N debe ser mayor a 0. Se ingresó =>%i<=\n", n);
        MPI_Finalize();
        return 1;    
    }
    int numThreads = atoi(argv[2]);
    if (numThreads <= 0){
        printf("El valor de numThreads debe ser mayor a 0. Se ingresó =>%i<=\n", numThreads);
        MPI_Finalize();
        return 1;    
    }
    int stripSize = n / numProcs;
    if (stripSize < blocksize){
        printf("Cambiaremos el blocksize al valor %i para una ejecución balanceada", stripSize);
        blocksize = stripSize;
    }

    if (n % blocksize != 0) {
        printf("El tamaño de la matriz (n x n) debe ser divisible por el tamaño del bloque\n");
        MPI_Finalize();
        return 2;
    }
    // ====================================VALIDACIONES===================================
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

    // Esto no cuenta como parte del tiempo? <===========================================================================================================================
    MPI_Scatter(A, bufferSizeStrip, MPI_DOUBLE, A, bufferSizeStrip, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Scatter(C, bufferSizeStrip, MPI_DOUBLE, C, bufferSizeStrip, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(B, size, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);

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
            // Obtener B transpuesta de todos los nodos
            MPI_Gather(BtransLoc, bufferSizeStrip, MPI_DOUBLE, BtransTot, bufferSizeStrip, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
            // Darles a todos el nuevo BT
            MPI_Bcast(BtransTot, n * n, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
            // Reducción a un valor del máximo, mínimo y sumas locales de A y B
            MPI_Reduce(&localMaxA, &maxA, 1, MPI_DOUBLE, MPI_MAX, MASTER, MPI_COMM_WORLD);
            MPI_Reduce(&localMinA, &minA, 1, MPI_DOUBLE, MPI_MIN, MASTER, MPI_COMM_WORLD);
            MPI_Reduce(&localSumaA, &sumaA, 1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
            MPI_Reduce(&localMaxB, &maxB, 1, MPI_DOUBLE, MPI_MAX, MASTER, MPI_COMM_WORLD);
            MPI_Reduce(&localMinB, &minB, 1, MPI_DOUBLE, MPI_MIN, MASTER, MPI_COMM_WORLD);
            MPI_Reduce(&localSumaB, &sumaB, 1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
            
            // Con los datos exactos, realizar el calculo del escalar
            if (rank == MASTER) {
                promedioA = sumaA / (n * n);
                promedioB = sumaB / (n * n);
                escalar = ((maxA * maxB) - (minA * minB)) / (promedioA * promedioB);
            }
            
            // Informar a todos los procesos del escalar
            MPI_Bcast(&escalar, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
        }

        multiplacionMatricesBloque(A, B, res1, blockSize, n, stripSize);
        multiplacionMatricesBloque(C, BtransTot, res2, blockSize, n, stripSize);
        multiplicarMatrizNumero(res1, escalar, stripSize, n);
        sumaMatrices(res1, res2, res1, n, stripSize);
    }
    
    // Obtener el resultado final a partir de lo parcial de todos los nodos
    MPI_Gather(res1, bufferSizeStrip, MPI_DOUBLE, R, bufferSizeStrip, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);

    
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
