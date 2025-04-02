#include "utils.h"
#include <stdint.h>

/*  Retorna tempo em milisegundos

    Forma de uso:

    double tempo;
    tempo = timestamp();
    <trecho de programa do qual se deseja medir tempo>
    tempo = timestamp() - tempo;
*/

double timestamp(void) {
  struct timespec tp;
  clock_gettime(CLOCK_MONOTONIC_RAW, &tp);
  return ((double)(tp.tv_sec * 1.0e3 + tp.tv_nsec * 1.0e-6));
}

real_t abs_value(real_t x) { return (x < 0) ? -x : x; }

int64_t ulp_distance(Double_t x, Double_t y) {
  return (abs_value((x.i - y.i)));
}
