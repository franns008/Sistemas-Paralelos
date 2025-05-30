#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <omp.h>
#include "mpi.h"

#define MASTER 0
#define tamBloque 128  // Tamaño de bloque para la multiplicación por bloques

void multiplicarMatricesBloquesColumnasBThread(double *A, double *B, double *C, int N, long inicio, long fin) {
    int indexII, indexJJ;
    
    #pragma omp for private(indexII, indexJJ) schedule(static)
    for (int i = inicio; i < fin; i += tamBloque) {
        for (int j = 0; j < N; j += tamBloque) {
            for (int k = 0; k < N; k += tamBloque) {
                for (int ii = i; ii < i + tamBloque && ii < fin; ii++) {
                    indexII = ii * N;
                    for (int jj = j; jj < j + tamBloque && jj < N; jj++) {
                        double aux = 0.0;
                        indexJJ = jj * N;
                        for (int kk = k; kk < k + tamBloque && kk < N; kk++) {
                            aux += A[indexII + kk] * B[indexJJ + kk];
                        }
                        C[indexII + jj] += aux;
                    }
                }
            }
        }
    }
}

int main(int argc, char *argv[]) {
    double *a, *b, *c, *d, *e, *f, *res1, *res2, *res3;
    int i, j, k, n, rank, numProcs, check = 1;
    MPI_Status status;
    double commTimes[4] = {0}, maxCommTimes[4] = {0}, minCommTimes[4] = {0}, commTime, totalTime;

    if ((argc != 2) || ((n = atoi(argv[1])) <= 0)) {
        printf("\nUsar: %s size \n  size: Dimension de la matriz y el vector\n", argv[0]);
        exit(1);
    }

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int stripSize = n / numProcs;

    // Asignación de memoria
    if (rank == MASTER) {
        a = (double *)malloc(sizeof(double) * n * n);
        c = (double *)malloc(sizeof(double) * n * n);
        f = (double *)malloc(sizeof(double) * n * n);
        d = (double *)malloc(sizeof(double) * n * n);
        e = (double *)malloc(sizeof(double) * n * n);
        b = (double *)malloc(sizeof(double) * n * n);
        res1 = (double *)malloc(sizeof(double) * n * n);
        res2 = (double *)malloc(sizeof(double) * n * n);
        res3 = (double *)malloc(sizeof(double) * n * n);
    } else {
        a = (double *)malloc(sizeof(double) * n * stripSize);
        c = (double *)malloc(sizeof(double) * n * stripSize);
        f = (double *)malloc(sizeof(double) * n * stripSize);
        b = (double *)malloc(sizeof(double) * n * n);
        d = (double *)malloc(sizeof(double) * n * n);
        e = (double *)malloc(sizeof(double) * n * n);
        res1 = (double *)calloc(n * stripSize, sizeof(double));
        res2 = (double *)calloc(n * stripSize, sizeof(double));
        res3 = (double *)calloc(n * stripSize, sizeof(double));
    }

    // Inicializar matrices en el master
    if (rank == MASTER) {
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                a[i * n + j] = 1;
                c[i * n + j] = 1;
                f[i * n + j] = 1;
                b[j * n + i] = 1;
                d[j * n + i] = 1;
                e[j * n + i] = 1;
                res1[i * n + j] = 0;
                res2[i * n + j] = 0;
                res3[i * n + j] = 0;
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    commTimes[0] = MPI_Wtime();
    // Distribuir datos
    MPI_Scatter(a, n * stripSize, MPI_DOUBLE, a, n * stripSize, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Scatter(c, n * stripSize, MPI_DOUBLE, c, n * stripSize, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Scatter(f, n * stripSize, MPI_DOUBLE, f, n * stripSize, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(b, n * n, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(d, n * n, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(e, n * n, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);

    commTimes[1] = MPI_Wtime();
    #pragma omp parallel 
    {
            // Inicializar matrices en los esclavos
        // Multiplicaciones + suma
        multiplicarMatricesBloquesColumnasBThread(a, b, res1, n, 0, stripSize);
        multiplicarMatricesBloquesColumnasBThread(c, d, res2, n, 0, stripSize);
        multiplicarMatricesBloquesColumnasBThread(f, e, res3, n, 0, stripSize);

        #pragma omp for private(i, j)
        for (i = 0; i < stripSize; i++) {
            for (j = 0; j < n; j++) {
                res1[i * n + j] += res2[i * n + j] + res3[i * n + j];
            }
        }
    }
    commTimes[2] = MPI_Wtime();
    // Recolectar resultados
    MPI_Gather(res1, n * stripSize, MPI_DOUBLE, res1, n * stripSize, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    
    // Sincronizar tiempos
    commTimes[3] = MPI_Wtime();
    MPI_Reduce(commTimes, minCommTimes, 4, MPI_DOUBLE, MPI_MIN, MASTER, MPI_COMM_WORLD);
    MPI_Reduce(commTimes, maxCommTimes, 4, MPI_DOUBLE, MPI_MAX, MASTER, MPI_COMM_WORLD);
 

    

    if (rank == MASTER) {
        int error_i = -1, error_j = -1;

        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                if (res1[i * n + j] != 3 * n) {
                    check = 0;
                    error_i = i;
                    error_j = j;
                    break;
                }
            }
            if (!check) break;
        }

        if (check) {
            printf("Multiplicacion de funcion resultado correcto\n");
        } else {
            printf("Multiplicacion incorrecta en res1[%d][%d] = %lf\n", error_i, error_j, res1[error_i * n + error_j]);
        }

        totalTime = maxCommTimes[3] - minCommTimes[0];
		commTime = (maxCommTimes[1] - minCommTimes[0]) + (maxCommTimes[3] - minCommTimes[2]);		

        printf("Multiplicacion (N=%d)\tTiempo total=%lf\tTiempo comunicacion=%lf\n", n, totalTime, commTime);
    }

    // Liberar memoria
    free(a); free(b); free(c); free(d); free(e); free(f);
    free(res1); free(res2); free(res3); 

    MPI_Finalize();
    
    return 0;
}
