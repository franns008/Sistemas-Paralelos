#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <sys/time.h>

#define T 10 // Número de hilos
long N = 500000000; // Tamaño del vector

double *vectorB;
double *vectorC;
double *vectorA;

// Función para medir el tiempo
double dwalltime() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec / 1000000.0;
}

// Inicializa un vector con un valor dado
void inicializarVector(double *vector, double val) {
    for (long i = 0; i < N; i++) {
        vector[i] = val;
    }
}

// Función que ejecutan los hilos para sumar los vectores
void *sumaVectorThread(void *arg) {
    long thread_id = *(long *)arg;
    free(arg);  // Liberamos la memoria asignada a thread_id

    long inicio = (N / T) * thread_id;
    long fin = (thread_id == T - 1) ? N : inicio + (N / T);

    for (long i = inicio; i < fin; i++) {
        vectorA[i] = vectorC[i] + vectorB[i];
    }
    printf("Hilo %d ha terminado su trabajo.\n", thread_id);
    pthread_exit(NULL);
}

// Función para imprimir el vector
void printVector(double *vector) {
    for (long i = 0; i < N; i++) {
        printf("%f ", vector[i]);
    }
    printf("\n");
}

long main() {
    pthread_t threads[T];

    // Asignar memoria dinámica
    vectorA = (double *)malloc(N * sizeof(double));
    vectorB = (double *)malloc(N * sizeof(double));
    vectorC = (double *)malloc(N * sizeof(double));

    if (!vectorA || !vectorB || !vectorC) {
        printf("Error al asignar memoria.\n");
        return 1;
    }

    // Inicializar vectores
    inicializarVector(vectorA, 0.0);
    inicializarVector(vectorB, 4.0);
    inicializarVector(vectorC, 2.0);

    // Crear los hilos
    for (long i = 0; i < T; i++) {
        long *thread_id = malloc(sizeof(long));
        *thread_id = i;
        pthread_create(&threads[i], NULL, sumaVectorThread, thread_id);
    }
    printf("Hilos creados. Esperando a que terminen...\n");
    double start_time = dwalltime();

    // Esperar a que todos los hilos terminen
    for (long i = 0; i < T; i++) {
        pthread_join(threads[i], NULL);
    }
    printf("Tiempo de ejecución con N = %d: %f segundos\n", N, dwalltime() - start_time);
    printf("Todos los hilos han terminado.\n");

    // Imprimir resultado
    //printVector(vectorA);

    // Liberar memoria
    free(vectorA);
    free(vectorB);
    free(vectorC);

    return 0;
}
