#include "bisection_func.h"
#include "transcendental.h"
#include <math.h>

double BisectionMethod_function(func_t *func, double a, double b, double *e, int *iterations)
{
  double c = (b + a) / 2;

  (*iterations)++;

  if (fabs(b - a) < 2 * *e || CountFunction(func, c) == 0)
  {
    return c;
  }
  else
  {
    if (CountFunction(func, c) * CountFunction(func, a) < 0.0)
    {
      return BisectionMethod_function(func, a, c, e, iterations);
    }
    else
    {
      return BisectionMethod_function(func, c, b, e, iterations);
    }
  }
}