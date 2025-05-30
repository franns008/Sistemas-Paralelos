#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include "mpi.h" 

#define MASTER 0

int tamBloque = 128;

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
                            aux += A[indexJJ + kk] * B[indexJJ + kk];
                        }
                    
                        C[indexII + jj] += aux; // Asignar el resultado a C
                        
                    }
                }
            }
        }
    }
}


int main(int argc, char *argv[]){
    
    MPI_Init(&argc,&argv);
    double *a, *b, *c, *d, *e, *f, *res1, *res2, *res3, *R;
    int i, j, k, n, rank, numProcs, check=1;
    MPI_Status status;
    double commTimes[4], maxCommTimes[4], minCommTimes[4], commTime, totalTime;
    
    if ((argc != 2) || ((n = atoi(argv[1])) <= 0) ) {
	    printf("\nUsar: %s size \n  size: Dimension de la matriz y el vector\n", argv[0]);
		exit(1);
	}
    
   
    MPI_Comm_size(MPI_COMM_WORLD,&numProcs);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    // calcular porcion de cada worker
    int stripSize = n / numProcs;
    printf("numProcs: %d, rank: %d, stripSize: %d\n", numProcs, rank, stripSize);
    // Reservar memoria
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
        R = (double *)malloc(sizeof(double) * n * n);
    } else {
        a = (double *)malloc(sizeof(double) * n * stripSize);
        c = (double *)malloc(sizeof(double) * n * stripSize);
        e = (double *)malloc(sizeof(double) * n * stripSize);
        b = (double *)malloc(sizeof(double) * n * n);
        d = (double *)malloc(sizeof(double) * n * n);
        f = (double *)malloc(sizeof(double) * n * n);
        res1 = (double *)calloc(n * stripSize, sizeof(double));
        res2 = (double *)calloc(n * stripSize, sizeof(double));
        res3 = (double *)calloc(n * stripSize, sizeof(double));
        R = (double *)malloc(sizeof(double) * n * stripSize);
    }
    // inicializar datos
    if (rank == 0) {
        for (i=0; i<n ; i++){
            for (j=0; j<n ; j++){
                a[i*n+j] = 1;
                c[i*n+j] = 1;
                e[i*n+j] = 1;
                res1[i*n+j] = 0;
                res2[i*n+j] = 0;
                res3[i*n+j] = 0;
            }
        }
        for (i=0; i<n ; i++){
            for (j=0; j<n ; j++){
                b[j*n+i] = 1;
                d[j*n+i] = 1;
                f[j*n+i] = 1;
            }
        }
    }else{
        for(i=0; i<n ; i++){
            for (j=0; j<n ; j++){
                res1[i*n+j] = 0;
                res2[i*n+j] = 0;
                res3[i*n+j] = 0;
            }
        }
    }
    printf("res1:%d res2:%d res3:%d\n",res1[1],res2[1],res3[1]);
    MPI_Barrier(MPI_COMM_WORLD);
    commTimes[0] = MPI_Wtime();
    // Broadcast the matrix and vector to all processes
    MPI_Scatter(a, n*stripSize, MPI_DOUBLE, a, n*stripSize, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Scatter(c, n*stripSize, MPI_DOUBLE, c, n*stripSize, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Scatter(e, n*stripSize, MPI_DOUBLE, e, n*stripSize, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(b, n*n, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(d, n*n, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(f, n*n, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);

    commTimes[1] = MPI_Wtime();
  
    

    multiplicarMatricesBloquesColumnasBThread(a, b, res1, n, 0, stripSize);
    multiplicarMatricesBloquesColumnasBThread(c, d, res2, n, 0, stripSize);
    multiplicarMatricesBloquesColumnasBThread(e, f, res3, n, 0, stripSize);
    
    /*
    Recolectar resultados parciales
    En este caso, cada proceso tiene su propia porciÃ³n de la matriz resultante
    y la matriz resultante completa se almacena en el proceso maestro
    MPI_Gather(res1, n*stripSize, MPI_DOUBLE, res1, n*stripSize, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Gather(res2, n*stripSize, MPI_DOUBLE, res2, n*stripSize, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Gather(res3, n*stripSize, MPI_DOUBLE, res3, n*stripSize, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);

    computar la suma de las matrices resultantes
    MPI_Scatter(res1, n*stripSize, MPI_DOUBLE, res1, n*stripSize, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Scatter(res2, n*stripSize, MPI_DOUBLE, res2, n*stripSize, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Scatter(res3, n*stripSize, MPI_DOUBLE, res3, n*stripSize, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    
    No uses MPI_Scatter despuÃ©s del Gather: eso destruye datos ya reunidos.
    
    Este bloque sobrescribe los datos que acabas de recolectar con MPI_Gather. Â¡Esto es incorrecto! 
    Lo que haces es:

    MPI_Gather: Traes los resultados parciales de todos los procesos al master.

    Luego haces MPI_Scatter: Sobrescribes esos resultados 
    (Â¡incluyendo en el master!) con porciones de lo que ya habÃ­as calculado 
    localmente antes de reunirlos.
    
    */
  
    printf("res1:%d res2:%d res3:%d\n",res1[1],res2[1],res3[1]);

    for(i=0; i<stripSize; i++) {
        for (j=0; j<n ;j++ ) {
            R[i*n + j] = res1[i*n+j] + res2[i*n+j] + res3[i*n+j];
        }
    }

    commTimes[2] = MPI_Wtime();
    

    MPI_Gather(R, n*stripSize, MPI_DOUBLE, R, n*stripSize, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    commTimes[3] = MPI_Wtime();

    MPI_Reduce(commTimes, minCommTimes, 4, MPI_DOUBLE, MPI_MIN, MASTER, MPI_COMM_WORLD);
    MPI_Reduce(commTimes, maxCommTimes, 4, MPI_DOUBLE, MPI_MAX, MASTER, MPI_COMM_WORLD);

   
    
    // Check results
    
    if (rank==MASTER) {
        int error_i, error_j;
		// Check results
        for(i = 0; i < n; i++) {
            for(j = 0; j < n; j++) {
                if (R[i*n + j] != 3 * n) {
                    check = 0;
                    error_i = i;
                    error_j = j;
                    printf("res1[%d][%d] = %lf\n", error_i, error_j, R[error_i*n+error_j]);
                    check = 1;
                }
            }
        }

		if(check){
			printf("Multiplicacion de funcion resultado correcto\n");
		}else{
			printf("Multiplicacion de funcion incorrecto, resultados obtenidos:\n");
            printf("res1[%d][%d] = %lf\n", error_i, error_j, res1[error_i*n+error_j]);
            
		}
		
		totalTime = maxCommTimes[3] - minCommTimes[0];
		commTime = (maxCommTimes[1] - minCommTimes[0]) + (maxCommTimes[3] - minCommTimes[2]);		

		printf("Multiplicacion de funcion (N=%d)\tTiempo total=%lf\tTiempo comunicacion=%lf\n",n,totalTime,commTime);
	}
   
    // Free memory
    free(a);
    free(b);
    free(c);
    free(d);
    free(e);
    free(f);
    free(res1);
    free(res2);
    free(res3);
    
    MPI_Finalize();    
    return 0;
}