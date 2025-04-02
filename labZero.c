#include "ZeroFuncao.h"
#include "utils.h"
#include <float.h>
#include <stdio.h>

int main() {
  int degree;
  Double_t a, b;
  Polinomio pol1;

  // Creates the first polynomial
  scanf("%d", &degree);
  pol1 = create_poly(degree);
  scanf("%lf %lf", &a.f, &b.f); // intervalo onde está uma das raizes.

  // #### START FAST #####
  printf("RAPIDO\n\n");

  // Setup bissec lento
  for (int i = 1; i <= 3; i++) {
    Double_t error;
    Double_t root;

    error.f = 0.0;
    root.f = 0.0;
    int n_it = 0;
    double time = 0.0;

    // Exec and run time
    time = timestamp();
    error.f = bisseccao(pol1, a, b, i, &n_it, calcPolinomio_rapido, &root);
    time = timestamp() - time;

    // bissec <raiz> <valor_crit-01> <it> <tempo>
    printf("bissec %.15e %.15e %d %.8e\n", root.f, error.f, n_it, time);
  }

  // Setup newton
  for (int i = 1; i <= 3; i++) {
    Double_t error;
    Double_t root;

    error.f = 0.0;
    root.f = 0.0;
    int n_it = 0;
    double time = 0.0;

    // Exec and run time
    time = timestamp();
    error.f = newtonRaphson(pol1, b, i, &n_it, calcPolinomio_rapido, &root);
    time = timestamp() - time;

    // newton <raiz> <valor_crit-01> <it> <tempo>
    printf("newton %.15e %.15e %d %.8e\n", root.f, error.f, n_it, time);
  }

  // #### END FAST #####
  printf("\nLENTO\n\n");
  // #### START SLOW #####

  // Setup bissec lento
  for (int i = 1; i <= 3; i++) {
    Double_t error;
    Double_t root;

    error.f = 0.0;
    root.f = 0.0;
    int n_it = 0;
    double time = 0.0;

    // Exec and run time
    time = timestamp();
    error.f = bisseccao(pol1, a, b, i, &n_it, calcPolinomio_lento, &root);
    time = timestamp() - time;

    // bissec <raiz> <valor_crit-01> <it> <tempo>
    printf("bissec %.15e %.15e %d %.8e\n", root.f, error.f, n_it, time);
  }

  // Setup newton
  for (int i = 1; i <= 3; i++) {
    Double_t error;
    Double_t root;

    error.f = 0.0;
    root.f = 0.0;
    int n_it = 0;
    double time = 0.0;

    // Exec and run time
    time = timestamp();
    error.f = newtonRaphson(pol1, b, i, &n_it, calcPolinomio_lento, &root);
    time = timestamp() - time;

    // newton <raiz> <valor_crit-01> <it> <tempo>
    printf("newton %.15e %.15e %d %.8e\n", root.f, error.f, n_it, time);
  }
  // #### END SLOW #####

  free_poly(&pol1);

  return 0;
}
