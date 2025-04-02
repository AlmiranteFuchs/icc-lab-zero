#ifndef __ZEROFUNCAO_H__
#define __ZEROFUNCAO_H__

#include "utils.h"
#include <float.h>

// Aproximação aceitável como valor zero
#define ZERO DBL_EPSILON

// Parâmetros para teste de convergência
#define MAXIT 500
#define EPS 1.0e-6
#define ULPS 2

typedef struct {
  Double_t *p;
  int grau;
} Polinomio;

Polinomio create_poly(int degree);
void free_poly(Polinomio *poly);

// Métodos
// Retornam valor do erro quando método finalizou. Este valor depende de
// tipoErro

real_t newtonRaphson(Polinomio p, Double_t x0, int criterioParada, int *it,
                     void (*calcFunc)(Polinomio, real_t, real_t *, real_t *),
                     Double_t *raiz);
real_t bisseccao(Polinomio p, Double_t a, Double_t b, int criterioParada,
                 int *it,
                 void (*calcFunc)(Polinomio, real_t, real_t *, real_t *),
                 Double_t *raiz);

// Cálculo de Polinômios
void calcPolinomio_rapido(Polinomio p, real_t x, real_t *px, real_t *dpx);
void calcPolinomio_lento(Polinomio p, real_t x, real_t *px, real_t *dpx);

#endif // __ZEROFUNCAO_H__
