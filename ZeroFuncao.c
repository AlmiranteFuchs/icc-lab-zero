#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "ZeroFuncao.h"
#include "utils.h"

// Creates the structure
Polinomio create_poly(int degree) {
  Polinomio poly;

  poly.grau = degree;
  poly.p = (Double_t *)malloc((poly.grau + 1) * sizeof(Double_t));

  if (poly.p == NULL) {
    printf("Memory alloc failed");
    exit(1);
  }

  // Reads the coefficients
  for (int i = poly.grau; i >= 0; i--) {
    scanf("%lf", &poly.p[i].f);
  }

  return poly;
}

// Clears the structure
void free_poly(Polinomio *poly) {
  if (poly->p != NULL) {
    free(poly->p);
    poly->p = NULL;
  }
}

// Retorna valor do erro quando método finalizou. Este valor depende de tipoErro
real_t newtonRaphson(Polinomio p, Double_t x0, int criterioParada, int *it,
                     void (*calcFunc)(Polinomio, real_t, real_t *, real_t *),
                     Double_t *raiz) {

  Double_t xk;
  Double_t xk_old;
  Double_t fxk;
  Double_t dpxk;
  Double_t L;

  xk.f = x0.f;
  xk_old.f = INFINITY;
  fxk.f = 0.0;
  dpxk.f = x0.f;
  L.f = 0.0;

  do {
    // Calc phi
    calcFunc(p, xk.f, &fxk.f, &dpxk.f);

    // Avoids div zero
    L.f = (abs_value(dpxk.f) > ZERO) ? dpxk.f : L.f;

    xk.f = xk.f - (fxk.f / L.f);

    *(raiz) = xk;

    // Circuit breaker
    if (criterioParada == 1) {
      if (abs_value(xk.f - xk_old.f) <= EPS) {
        return abs_value(xk.f - xk_old.f);
      }
    }

    if (criterioParada == 2) {
      if ((abs_value(fxk.f)) <= ZERO) {
        return abs_value(fxk.f);
      }
    }

    if (criterioParada == 3) {
      uint64_t ulps = ulp_distance(xk_old, xk);

      if (ulps <= 2) {
        return xk.f;
      }
    }

    xk_old = xk;
    *(it) = *(it) + 1;

  } while (*it < MAXIT);
}

// Retorna valor do erro quando método finalizou. Este valor depende de tipoErro
real_t bisseccao(Polinomio p, Double_t a, Double_t b, int criterioParada,
                 int *it,
                 void (*calcFunc)(Polinomio, real_t, real_t *, real_t *),
                 Double_t *raiz) {
  Double_t xl = a;
  Double_t xu = b;
  Double_t xm;
  Double_t xm_old;
  Double_t fxm;
  Double_t fxl;
  Double_t dpxm;
  Double_t dpxl;

  xm.f = 0.0;
  xm_old.f = INFINITY;
  fxm.f = 0.0;
  fxl.f = 0.0;
  dpxm.f = 0.0;
  dpxl.f = 0.0;

  do {
    xm.f = ((xl.f + xu.f) / 2);

    calcFunc(p, xl.f, &fxl.f, &dpxl.f);
    calcFunc(p, xm.f, &fxm.f, &dpxm.f);

    if (fxl.f * fxm.f < 0) {
      // The root is between this interval
      xu = xm;

    } else {
      // The root is not between this interval
      xl = xm;
    }

    *(raiz) = xm;

    // Circuit breaker
    if (criterioParada == 1) {
      if (abs_value(xm.f - xm_old.f) <= EPS) {
        return abs_value(xm.f - xm_old.f);
      }
    }

    if (criterioParada == 2) {
      if ((abs_value(fxm.f)) <= ZERO) {
        return abs_value(fxm.f);
      }
    }

    if (criterioParada == 3) {
      int ulps = ulp_distance(xm_old, xm);

      if (ulps <= 2) {
        return xm.f;
      }
    }

    xm_old = xm;
    *(it) = *(it) + 1;

  } while (*it < MAXIT);
}

void calcPolinomio_rapido(Polinomio p, real_t x, real_t *px, real_t *dpx) {
  real_t b = 0;
  real_t c = 0;

  for (int i = p.grau; i > 0; i--) {
    b = b * x + p.p[i].f;
    c = c * x + b;
  }

  b = b * x + p.p[0].f;
  *px = b;
  *dpx = c;
}

void calcPolinomio_lento(Polinomio p, real_t x, real_t *px, real_t *dpx) {
  *px = p.p[0].f;
  *dpx = 0.0;

  for (int i = 1; i <= p.grau; i++) {
    *px += (p.p[i].f) * (pow(x, i));
    *dpx += i * (p.p[i].f * pow(x, i - 1));
  }
}