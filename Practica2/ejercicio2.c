#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

#define T 4 // Número de hilos

double *A;
double *B;
double *C1, *C2;
int N = 16; // Tamaño de la matriz
int bs = 4; // Tamaño del bloque

pthread_mutex_t *mutexVec;

void initvalmat(double *mat, int n, double val, int transpose){
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (transpose == 0)
                mat[i * n + j] = val;
            else
                mat[j * n + i] = val;
        }
    }
}

void blkmul(double *ablk, double *bblk, double *cblk, int n, int bs) {
    for (int i = 0; i < bs; i++) {
        for (int j = 0; j < bs; j++) {
            for (int k = 0; k < bs; k++) {
                cblk[i * n + j] += ablk[i * n + k] * bblk[k * n + j];
            }
        }
    }
}

void blkmulThread(double *ablk, double *bblk, double *cblk, int n, int bs, int row_offset, int col_offset) {
    for (int i = 0; i < bs; i++) {
        for (int j = 0; j < bs; j++) {
            for (int k = 0; k < bs; k++) {
                int index = (row_offset + i) * n + (col_offset + j);
                //pthread_mutex_lock(&mutexVec[index]);
                cblk[i * n + j] += ablk[i * n + k] * bblk[k * n + j];
                //pthread_mutex_unlock(&mutexVec[index]);
            }
        }
    }
}

void matmulblksThreads(double *a, double *b, double *c, int n, int bs, long block_start, long block_end) {
    for (int bi = block_start; bi < block_end; bi++) {
        int i = (bi * bs);
        for (int j = 0; j < n; j += bs) {
            for (int k = 0; k < n; k += bs) {
                blkmulThread(&a[i * n + k], &b[k * n + j], &c[i * n + j], n, bs, i, j);
            }
        }
    }
}

void matmulblks(double *a, double *b, double *c, int n, int bs) {
    for (int i = 0; i < n; i += bs) {
        for (int j = 0; j < n; j += bs) {
            for (int k = 0; k < n; k += bs) {
                blkmul(&a[i * n + k], &b[k * n + j], &c[i * n + j], n, bs);
            }
        }
    }
}

void *multiVectorThread(void *arg) {
    long thread_id = *(long *)arg;
    free(arg);

    long total_blocks = N / bs;
    long blocks_per_thread = total_blocks / T;
    long block_start = thread_id * blocks_per_thread;
    long block_end = (thread_id == T - 1) ? total_blocks : block_start + blocks_per_thread;

    matmulblksThreads(A, B, C1, N, bs, block_start, block_end);

    pthread_exit(NULL);
}

void printMat(double *mat, int N) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%5.1f ", mat[i * N + j]);
        }
        printf("\n");
    }
    printf("\n");
}

void validar(double *mat1, double *mat2){
    for (int i = 0; i < N * N; i++) {
        if (mat1[i] != mat2[i]) {
            printf("Error en la validacion en posición %d\n", i);
            return;
        }
    }
    printf("El calculo de la ecuacion es correcto\n");
}

int main() {
    pthread_t threads[T];

    // Asignar memoria
    A = (double *)malloc(N * N * sizeof(double));
    B = (double *)malloc(N * N * sizeof(double));
    C1 = (double *)malloc(N * N * sizeof(double));
    C2 = (double *)malloc(N * N * sizeof(double));
    mutexVec = (pthread_mutex_t *)malloc(N * N * sizeof(pthread_mutex_t));

    if (!A || !B || !C1 || !C2 || !mutexVec) {
        printf("Error al asignar memoria.\n");
        return 1;
    }

    // Inicializar matrices
    initvalmat(C1, N, 0.0, 0);
    initvalmat(C2, N, 0.0, 0);
    initvalmat(A, N, 3.0, 0);
    initvalmat(B, N, 1.0, 1); // Transpuesta

    // Inicializar mutexes
    for (int i = 0; i < N * N; i++) {
        pthread_mutex_init(&mutexVec[i], NULL);
    }

    // Crear hilos
    for (long i = 0; i < T; i++) {
        long *thread_id = malloc(sizeof(long));
        *thread_id = i;
        pthread_create(&threads[i], NULL, multiVectorThread, thread_id);
    }

    // Esperar a que terminen
    for (int i = 0; i < T; i++) {
        pthread_join(threads[i], NULL);
    }

    printf("Resultado con hilos:\n");
    printMat(C1, N);

    matmulblks(A, B, C2, N, bs);
    printf("Resultado secuencial:\n");
    printMat(C2, N);

    validar(C1, C2);

    // Liberar recursos
    for (int i = 0; i < N * N; i++) {
        pthread_mutex_destroy(&mutexVec[i]);
    }
    free(mutexVec);
    free(A); free(B); free(C1); free(C2);
    return 0;
}
