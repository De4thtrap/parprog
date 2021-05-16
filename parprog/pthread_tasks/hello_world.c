#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>

#define DEFAULT_N 4

void *hello(void *param);


int main(int argc, char *argv[]) {

    // thread count
    long N = DEFAULT_N;
    if (argc >= 2) {
        N = strtol(argv[1], NULL, 10);
    }

    pthread_t *threads = calloc(N, sizeof(pthread_t));

    // initializing threads
    for (int i = 0; i < N; i++) {
        pthread_create(&threads[i], NULL, hello, NULL);
    }

    // waiting for finishing threads
    for (int i = 0; i < N; i++) {
        pthread_join(threads[i], NULL);
    }

    free(threads);
}

// function for running by threads
void *hello(void *param) {
    printf("Hello, world! pid = %ld\n", pthread_self());
    pthread_exit(0);
}

