#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include "mpi.h" 
#define EPSILON 1e-6
#define MASTER 0

int iguales(double a, double b) {
        return fabs(a - b) < EPSILON;
}


int main(int argc, char *argv[]){
    double *v;
    
    int check=1;
    
    int n, rank, numProcs, i, j, k;
    double minLocal = 99999999;
    double maxLocal = -1;
    double promedioLoc = 0;
    double minGlobal, maxGlobal, promedioGlob = 0;

    double intervaloMin = 1.5;
    double intervaloMax = 10000;

    MPI_Status status;
    double commTimes[4], maxCommTimes[4], minCommTimes[4], commTime, totalTime;

    if (argc != 2) {
        printf("Usa: %s <size>\n", argv[0]);
        exit(1);
    }
    n = atoi(argv[1]);
    if (n <= 0) {
        printf("El tamaño de la matriz debe ser mayor que 0.\n");
        exit(1);
    }

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    if (n % numProcs != 0) {
        printf("El tamaño de la matriz debe ser múltiplo del número de procesos.\n");
        exit(1);
    }

    // calcular porción de cada worker
    int stripSize = n / numProcs;
    
    // Reservar memoria
    if (rank == MASTER) {
        v = (double*) malloc(sizeof(double) * n); // El maestro necesita todo
    } else {
        v = (double*) malloc(sizeof(double) * stripSize); // Los demás solo su porción
    }
    // inicializar datos
    if(rank == MASTER) {
        srand(time(NULL));
        for (i = 0; i < n; i++) {
            double r = (double)rand() / RAND_MAX;
            v[i] = r * (intervaloMax - intervaloMin) + intervaloMin;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    commTimes[0] = MPI_Wtime();
    MPI_Scatter(v, stripSize, MPI_DOUBLE, v, stripSize, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    commTimes[1] = MPI_Wtime();

    for (i = 0; i < stripSize; i++) {
        if(v[i] < minLocal){
            minLocal = v[i];
        }
        if(v[i] > maxLocal){
           maxLocal = v[i];
        }
        promedioLoc += v[i];
    }

    commTimes[2] = MPI_Wtime();

    MPI_Reduce(&minLocal, &minGlobal, 1, MPI_DOUBLE, MPI_MIN, MASTER, MPI_COMM_WORLD);
    MPI_Reduce(&maxLocal, &maxGlobal, 1, MPI_DOUBLE, MPI_MAX, MASTER, MPI_COMM_WORLD);
    MPI_Reduce(&promedioLoc, &promedioGlob, 1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);

    commTimes[3] = MPI_Wtime();

    MPI_Reduce(commTimes, minCommTimes, 4, MPI_DOUBLE, MPI_MIN, MASTER, MPI_COMM_WORLD);
    MPI_Reduce(commTimes, maxCommTimes, 4, MPI_DOUBLE, MPI_MAX, MASTER, MPI_COMM_WORLD);

    MPI_Finalize();
    // Check results
    if (rank == MASTER) {
        promedioGlob = promedioGlob / n;
        printf("Minimo: %lf, Maximo: %lf, Promedio: %lf\n", minGlobal, maxGlobal, promedioGlob);
        totalTime = maxCommTimes[3] - minCommTimes[0];
        commTime = (maxCommTimes[1] - minCommTimes[0]) + (maxCommTimes[3] - minCommTimes[2]);		
        printf("Tiempo total=%lf\tTiempo comunicacion=%lf\n", totalTime, commTime);
        
        double checkMin = 99999999;
        double checkMax = -1;
        double checkProm = 0;
        for(i=0;i<n;i++){
            if(v[i] < checkMin){
                checkMin = v[i];
            }
            if(v[i] > checkMax){
                checkMax = v[i];
            }
            checkProm += v[i];
        }       
        checkProm /= n;


        // comparo los resultados asi porque con punto flotante no se puede comparar directamente
        if(iguales(checkMin, minGlobal) && iguales(checkMax, maxGlobal) && iguales(checkProm, promedioGlob)) {
            printf("Los resultados son correctos\n");
        } else {
            printf("Los resultados son erroneos\n");
            printf("Minimo: %lf, Maximo: %lf, Promedio: %lf\n", checkMin, checkMax, checkProm);
        }
    }
    
    
    // Broadcast the vector to all processes
    free(v);

    

}
