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
    double timetick[2],totalTime,commTime;
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
    int size = n * n;
    int bufferStripSize = n * stripSize;
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
        BtransTot = (double *)malloc(sizeof(double) * bufferStripSize);
        res1 = (double *)calloc(bufferStripSize, sizeof(double));
        res2 = (double *)calloc(bufferStripSize, sizeof(double));
    }
    R = (double *)malloc(sizeof(double) * size);

    if (rank == MASTER) {
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                A[i * n + j] = i*j;
                C[i * n + j] = i*j;
            }
        }
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                B[i * n + j] = 5;
            }
        }
      
    }
    // ===================== PRIMERA ZONA =============================
    if (rank == MASTER) {
        timetick[0] = dwalltime();
    }
    MPI_Scatter(A, bufferStripSize, MPI_DOUBLE, A, bufferStripSize, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Scatter(C, bufferStripSize, MPI_DOUBLE, C, bufferStripSize, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(B, size, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    if (rank == MASTER) {
        commTime = dwalltime() - timetick[0];
    
    }

    MPI_Barrier(MPI_COMM_WORLD);

     // ===================== FIN PRIMERA ZONA =============================
    calcularMaximoMinimoPromedio(A, n, stripSize, &localMax[0], &localMin[0], &localSuma[0]);
    calcularMaximoMinimoPromedio(B, n, stripSize, &localMax[1], &localMin[1], &localSuma[1]);


    if(rank == MASTER) {
        timetick[1] = dwalltime();
    }
    MPI_Reduce(&localMax, &max, 2, MPI_DOUBLE, MPI_MAX, MASTER, MPI_COMM_WORLD);
    MPI_Reduce(&localMin, &min, 2, MPI_DOUBLE, MPI_MIN, MASTER, MPI_COMM_WORLD);
    MPI_Reduce(&localSuma, &suma, 2, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD); 

    if(rank == MASTER) {
        commTime += (dwalltime() - timetick[1]);
    }
    
     // ===================== SEGUNDA ZONA =============================
    double promedioA, promedioB, escalar;
    
    if (rank == MASTER) {
        printf("El maximo de A es ==> %f \n", max[0]);
        printf("El maximo de B es ==> %f \n", max[1]);
        printf("El minimo de A es ==> %f \n", min[0]);
        printf("El minimo de B es ==> %f \n", min[1]);
        printf("La suma total de A es ==> %f \n", suma[0]);
        printf("La suma total de B es ==> %f \n", suma[1]);
        promedioA = suma[0] / (size);
        promedioB = suma[1] / (size);
        escalar = ((max[0] * max[1]) - (min[0] * min[1])) / (promedioA * promedioB);
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
