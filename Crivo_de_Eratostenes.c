#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void sieve_initial(int limit, int **primes, int *count) {
    char *is_prime = malloc((limit + 1) * sizeof(char));
    for (int i = 0; i <= limit; i++) is_prime[i] = 1;
    is_prime[0] = is_prime[1] = 0;

    for (int p = 2; p * p <= limit; p++) {
        if (is_prime[p]) {
            for (int i = p * p; i <= limit; i += p)
                is_prime[i] = 0;
        }
    }

    int c = 0;
    for (int i = 2; i <= limit; i++) if (is_prime[i]) c++;

    *primes = malloc(c * sizeof(int));
    int idx = 0;
    for (int i = 2; i <= limit; i++)
        if (is_prime[i]) (*primes)[idx++] = i;

    *count = c;
    free(is_prime);
}

void sieve_block(int start, int end, int *primes, int prime_count, int **result, int *res_count) {
    char *is_prime = malloc((end - start + 1) * sizeof(char));
    for (int i = 0; i <= end - start; i++) is_prime[i] = 1;

    for (int i = 0; i < prime_count; i++) {
        int p = primes[i];
        int first = ((start + p - 1)/p) * p;
        if (first < p*p) first = p*p;
        for (int j = first; j <= end; j += p)
            is_prime[j - start] = 0;
    }

    int count = 0;
    for (int i = 0; i <= end - start; i++) if (is_prime[i]) count++;

    *result = malloc(count * sizeof(int));
    int idx = 0;
    for (int i = 0; i <= end - start; i++)
        if (is_prime[i]) (*result)[idx++] = start + i;

    *res_count = count;
    free(is_prime);
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc != 2) {
        if (rank == 0) printf("Uso: %s N\n", argv[0]);
        MPI_Finalize();
        return 1;
    }

    int N = atoi(argv[1]);
    double t0, t1;

    int sqrtN = (int)sqrt(N);
    int *primes_sqrt = NULL;
    int prime_count = 0;

    if (rank == 0) {
        sieve_initial(sqrtN, &primes_sqrt, &prime_count);
        t0 = MPI_Wtime();
    }

    MPI_Bcast(&prime_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (rank != 0) primes_sqrt = malloc(prime_count * sizeof(int));
    MPI_Bcast(primes_sqrt, prime_count, MPI_INT, 0, MPI_COMM_WORLD);

    int num_slaves = size - 1;
    int start = 2 + (rank-1)*(N-1)/num_slaves;
    int end   = 2 + rank*(N-1)/num_slaves - 1;
    if (rank == size-1) end = N;

    int *local_primes = NULL;
    int local_count = 0;
    if (rank != 0) {
        sieve_block(start, end, primes_sqrt, prime_count, &local_primes, &local_count);
        MPI_Send(&local_count, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(local_primes, local_count, MPI_INT, 0, 1, MPI_COMM_WORLD);
        free(local_primes);
    } else {
        FILE *fout = fopen("primos.txt", "w");
        for (int i = 0; i < prime_count; i++) fprintf(fout, "%d\n", primes_sqrt[i]);
        for (int s = 1; s < size; s++) {
            int count;
            MPI_Recv(&count, 1, MPI_INT, s, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            int *arr = malloc(count * sizeof(int));
            MPI_Recv(arr, count, MPI_INT, s, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (int i = 0; i < count; i++) fprintf(fout, "%d\n", arr[i]);
            free(arr);
        }
        fclose(fout);
        t1 = MPI_Wtime();
        printf("Tempo de execução: %f segundos\n", t1-t0);
    }

    free(primes_sqrt);
    MPI_Finalize();
    return 0;
}
