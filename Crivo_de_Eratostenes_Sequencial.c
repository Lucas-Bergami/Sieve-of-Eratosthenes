#include <stdio.h>
#include <stdlib.h>
#include <math.h>

char* init_sieve(int N) {
  char *is_prime = malloc((N + 1) * sizeof(char));
  for (int i = 0; i <= N; i++) is_prime[i] = 1;
  is_prime[0] = is_prime[1] = 0;
  return is_prime;
}

/* Marca os múltiplos de cada primo */
void mark_multiples(int N, char *is_prime) {
  for (int p = 2; p * p <= N; p++) {
    if (is_prime[p]) {
      for (int i = p * p; i <= N; i += p)
        is_prime[i] = 0;
    }
  }
}

int* collect_primes(int N, char *is_prime, int *count) {
  int c = 0;
  for (int i = 2; i <= N; i++) if (is_prime[i]) c++;

  int *primes = malloc(c * sizeof(int));
  int idx = 0;
  for (int i = 2; i <= N; i++)
    if (is_prime[i]) primes[idx++] = i;

  *count = c;
  return primes;
}

int main(int argc, char **argv) {
  if (argc != 2) {
    printf("Uso: %s N\n", argv[0]);
    return 1;
  }

  int N = atoi(argv[1]);
  int count = 0;

  char *is_prime = init_sieve(N);

  mark_multiples(N, is_prime);

  int *primes = collect_primes(N, is_prime, &count);

  FILE *fout = fopen("primos.txt", "w");
  if (!fout) {
    perror("Erro ao abrir primos.txt");
    free(primes);
    free(is_prime);
    return 1;
  }

  for (int i = 0; i < count; i++){
    fprintf(fout, "%d", primes[i]);
    fprintf(fout, ",");
  }

  fclose(fout);
  free(primes);
  free(is_prime);

  printf("Encontrados %d primos até %d. Resultado em primos.txt\n", count, N);

  return 0;
}
