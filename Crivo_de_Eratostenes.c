#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <mpi.h>

typedef struct {
    uint64_t *primes;
    size_t count;
} PrimeList;

#define GET_BIT(arr, i)   (((arr)[(i) >> 3] >> ((unsigned)((i) & 7))) & 1U)
#define CLEAR_BIT(arr, i) ((arr)[(i) >> 3] &= (unsigned char)~(1U << ((unsigned)((i) & 7))))

PrimeList sieve_initial(uint64_t limit) {
    PrimeList result = {NULL, 0};
    if (limit < 2) return result;

    uint64_t num_odds = (limit - 1) / 2; // somente ímpares > 2
    size_t array_size = (num_odds + 8) / 8;
    unsigned char *is_prime = malloc(array_size);
    if (!is_prime) { perror("malloc"); exit(EXIT_FAILURE); }
    memset(is_prime, 0xFF, array_size);

    uint64_t sqrt_limit = (uint64_t)sqrt((long double)limit);
    for (uint64_t p = 3; p <= sqrt_limit; p += 2) {
        uint64_t idx_p = (p - 3) / 2;
        if (GET_BIT(is_prime, idx_p)) {
            uint64_t start = (p * p - 3) / 2;
            for (uint64_t i = start; i < num_odds; i += p)
                CLEAR_BIT(is_prime, i);
        }
    }

    size_t count = 1;
    for (size_t i = 0; i < array_size; i++)
        count += __builtin_popcount((unsigned)is_prime[i]);

    if (count > SIZE_MAX / sizeof(uint64_t)) {
        fprintf(stderr, "Overflow na alocação de primos\n"); free(is_prime); exit(EXIT_FAILURE);
    }

    result.primes = malloc(count * sizeof(uint64_t));
    if (!result.primes) { perror("malloc"); free(is_prime); exit(EXIT_FAILURE); }

    size_t idx = 0;
    result.primes[idx++] = 2;
    for (uint64_t i = 0; i < num_odds; i++)
        if (GET_BIT(is_prime, i))
            result.primes[idx++] = 2 * i + 3;

    result.count = count;
    free(is_prime);
    return result;
}

void sieve_block(uint64_t start, uint64_t end, uint64_t *base_primes, size_t prime_count,
                 uint64_t **result, size_t *res_count) {

    size_t size = end - start + 1;
    unsigned char *is_prime = malloc((size + 7) / 8);
    if (!is_prime) { perror("malloc"); exit(EXIT_FAILURE); }
    memset(is_prime, 0xFF, (size + 7) / 8);

    for (size_t i = 0; i < prime_count; i++) {
        uint64_t p = base_primes[i];
        if (p < 2) continue; // proteção
        uint64_t first = ((start + p - 1) / p) * p;
        if (first < p * p) first = p * p;
        for (uint64_t j = first; j <= end; j += p)
            CLEAR_BIT(is_prime, j - start);
    }

    size_t count = 0;
    for (uint64_t i = 0; i < size; i++)
        if (GET_BIT(is_prime, i)) count++;

    *result = malloc(count * sizeof(uint64_t));
    if (!*result) { perror("malloc"); free(is_prime); exit(EXIT_FAILURE); }

    size_t idx = 0;
    for (uint64_t i = 0; i < size; i++)
        if (GET_BIT(is_prime, i))
            (*result)[idx++] = start + i;

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

    uint64_t N = strtoull(argv[1], NULL, 10);
    double t0, t1;

    uint64_t sqrtN = (uint64_t)sqrt((long double)N);
    uint64_t *primes_sqrt = NULL;
    size_t prime_count = 0;

    if (rank == 0) {
        PrimeList base = sieve_initial(sqrtN);
        primes_sqrt = base.primes;
        prime_count = base.count;
        t0 = MPI_Wtime();
    }

    MPI_Bcast(&prime_count, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
    if (rank != 0) primes_sqrt = malloc(prime_count * sizeof(uint64_t));
    MPI_Bcast(primes_sqrt, prime_count, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);

    int num_slaves = size - 1;
    uint64_t start = 2 + (rank - 1) * (N - 1) / num_slaves;
    uint64_t end   = 2 + rank * (N - 1) / num_slaves - 1;
    if (rank == size - 1) end = N;

    uint64_t *local_primes = NULL;
    size_t local_count = 0;

    if (rank != 0) {
        sieve_block(start, end, primes_sqrt, prime_count, &local_primes, &local_count);
        MPI_Send(&local_count, 1, MPI_UNSIGNED_LONG_LONG, 0, 0, MPI_COMM_WORLD);
        MPI_Send(local_primes, local_count, MPI_UNSIGNED_LONG_LONG, 0, 1, MPI_COMM_WORLD);
        free(local_primes);
    } else {
        FILE *fout = fopen("primos.txt", "w");
        for (size_t i = 0; i < prime_count; i++)
            fprintf(fout, "%llu\n", (unsigned long long)primes_sqrt[i]);

        for (int s = 1; s < size; s++) {
            size_t count;
            MPI_Recv(&count, 1, MPI_UNSIGNED_LONG_LONG, s, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            uint64_t *arr = malloc(count * sizeof(uint64_t));
            MPI_Recv(arr, count, MPI_UNSIGNED_LONG_LONG, s, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (size_t i = 0; i < count; i++)
                fprintf(fout, "%llu\n", (unsigned long long)arr[i]);
            free(arr);
        }

        fclose(fout);
        t1 = MPI_Wtime();
        printf("Tempo de execução: %f segundos\n", t1 - t0);
    }

    free(primes_sqrt);
    MPI_Finalize();
    return 0;
}
