#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <sys/time.h>

#define T 2
int N=2;
double promedioT =0;
int minimoT = 999999999;
int maximoT = -1;
int cantHilosTerminaron = 0;
int *vectorA;
pthread_mutex_t mutexMin,mutexMax,mutexProm; // Declaración del mutex

double dwalltime() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec / 1000000.0;
}

void inicializarVector(int *vector) {
    for (long i = 0; i < N; i++) {
        vector[i] = rand() % 1000; // Inicializa con valores aleatorios
    }
}

void *recorrerVector(void *args) {
    long thread_id = *(long *)args;
    free(args);  // Liberamos la memoria asignada a thread_id

    int minimo= 999999999;
    int maximo= -1;
    int promedio=0;

    long inicio = (N / T) * thread_id;
    long fin = (thread_id == T - 1) ? N : inicio + (N / T);

    int suma = 0;
    for (long i = inicio; i < fin; i++) {
        suma += vectorA[i];
        if (vectorA[i] < minimo) {
            minimo = vectorA[i];
        }
        if (vectorA[i] > maximo) {
            maximo = vectorA[i];
        }
    }
    pthread_mutex_lock(&mutexMin); // Bloqueamos el mutex
    if (minimo < minimoT) {
        minimoT = minimo;
    }
    pthread_mutex_unlock(&mutexMin); // Desbloqueamos el mutex
    pthread_mutex_lock(&mutexMax); // Bloqueamos el mutex
    if (maximo > maximoT) {
        maximoT = maximo;
    }
    pthread_mutex_unlock(&mutexMax); // Desbloqueamos el mutex
    pthread_mutex_lock(&mutexProm); // Bloqueamos el mutex
    cantHilosTerminaron++;
    promedioT += suma;
    if (cantHilosTerminaron == T) {
        promedioT /= N; // Calculamos el promedio total
    }
    pthread_mutex_unlock(&mutexProm);
     // Desbloqueamos el mutex
    printf("Hilo %ld ha terminado su trabajo. Suma parcial: %d\n", thread_id, suma);
    pthread_exit(NULL);
}

void validarVector(int *vector) {
    int suma = 0;
    int minimo = 999999999;
    int maximo = -1;
    double promedio = 0.0;  
    for (long i = 0; i < N; i++) {
        suma += vector[i];
        if (vector[i] < minimo) {
            minimo = vector[i];
        }
        if (vector[i] > maximo) {
            maximo = vector[i];
        }
    }
    promedio = (double)suma / N;
    if (promedio != promedioT) {
        printf("Error: el promedio no coincide. Promedio sin hilos: %f, Promedio con hilos: %f\n", promedio, promedioT);
    } else {
        printf("El promedio coincide. Promedio: %f\n", promedio);
    }
    if (minimo != minimoT) {
        printf("Error: el mínimo no coincide. Mínimo sin hilos: %d, Mínimo con hilos: %d\n", minimo, minimoT);
    } else {
        printf("El mínimo coincide. Mínimo: %d\n", minimo);
    }
    if (maximo != maximoT) {
        printf("Error: el máximo no coincide. Máximo sin hilos: %d, Máximo con hilos: %d\n", maximo, maximoT);
    } else {
        printf("El máximo coincide. Máximo: %d\n", maximo);
    }
}

void main() {
    pthread_t threads[T];
    vectorA = (int *)malloc(N * sizeof(int));
    inicializarVector(vectorA);

    for (long i = 0; i < T; i++) {
        long *thread_id = malloc(sizeof(long));
        *thread_id = i;
        pthread_create(&threads[i], NULL, recorrerVector, thread_id);
    }

    for (long i = 0; i < T; i++) {
        pthread_join(threads[i], NULL);
    }

    validarVector(vectorA);

    free(vectorA); // Liberar la memoria asignada al vector
}