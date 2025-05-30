#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <sys/time.h>


#define T 10 // Número de hilos
long N = 50000; // Tamaño del vector
double ValorBuscado = 500; // Valor a buscar en el vector
double *vectorA;
int total = 0;
pthread_mutex_t mutex; // Declaración del mutex

double dwalltime() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec / 1000000.0;
}

void inicializarVector(double *vector) {
    for (long i = 0; i < N; i++) {
        vector[i] = rand() % 1000; // Inicializa con valores aleatorios
    }
}

void *encontrarValor(void *args){
    long thread_id = *(long *)args;
    free(args);  // Liberamos la memoria asignada a thread_id

    long inicio = (N / T) * thread_id;
    long fin = (thread_id == T - 1) ? N : inicio + (N / T);

    int cant = 0;
    for (long i =  inicio;i < fin; i++){
        if (vectorA[i] == ValorBuscado){
            cant++;
        }
    }
    pthread_mutex_lock(&mutex); // Bloqueamos el mutex
    total += cant; // Acumulamos el resultado
    pthread_mutex_unlock(&mutex); // Desbloqueamos el mutex
    printf("Hilo %ld ha terminado su trabajo. Encontró %d veces el valor %f.\n", thread_id, cant, ValorBuscado);
    
}
void validarDatos(double *vector){
    int cantSinThread = 0;
    for (long i = 0; i < N; i++) {
        if (vector[i] == ValorBuscado) {
            cantSinThread++;                
        }
    }
    if (cantSinThread != total) {
        printf("Error: el conteo no coincide. Se encontraron %d veces el valor %f sin hilos y %d veces con hilos.\n", cantSinThread, ValorBuscado, total);
    } else {
        printf("El conteo coincide. Se encontraron %d veces el valor %f.\n", cantSinThread, ValorBuscado);
    }
}
void main(){
    pthread_t threads[T];
    vectorA = (double *)malloc(N * sizeof(double));
    inicializarVector(vectorA);
    pthread_mutex_init(&mutex, NULL); // Inicializamos el mutex

    for (long i = 0; i < T; i++) {
        long *thread_id = malloc(sizeof(long));
        *thread_id = i;
        pthread_create(&threads[i], NULL,encontrarValor, thread_id);
    }
    
    for (long i = 0; i < T; i++) {
        pthread_join(threads[i], NULL);
    }
    printf("El valor %f fue encontrado %d veces en el vector.\n", ValorBuscado, total);

    pthread_mutex_destroy(&mutex); // Destruimos el mutex
    validarDatos(vectorA); // Validamos los datos
    // Liberar memoria
    free(vectorA);

}
