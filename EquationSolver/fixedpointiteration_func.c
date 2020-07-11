#include "fixedpointiteration_func.h"
#include "transcendental.h"
#include <stdio.h>
#include <math.h>
#define CHECK_INTERVAL
#define FIRST_ANSWERS

static int _checkIntervalFunc(func_t *func, double a, double b)
{
  if (CountFunction(func, a) * CountFunction(func, b) < 0.0)
  {
    return 1;
  }
  else
  {
    return 0;
  }
}

static int _checkIntervalDerivativeFunc(func_t *der, double a, double b)
{
  double delta = 0.01;

  if (CountFunction(der, a) * CountFunction(der, b) < 0.0)
  {
    return 0;
  }
  else
  {
    while (a + delta < b)
    {
      if (CountFunction(der, a) * CountFunction(der, a + delta) < 0.0)
      {
        return 0;
      }
      a += delta;
    }
    return 1;
  }
}

double FindMaxAbsFunctionBorders(func_t *func, double a, double b)
{
  double aX = fabs(CountFunction(func, a));
  double bX = fabs(CountFunction(func, b));

  if (aX > bX)
  {
    return aX;
  }
  else
  {
    return bX;
  }
}

double FindMinAbsFunctionBorders(func_t *func, double a, double b)
{
  double aX = fabs(CountFunction(func, a));
  double bX = fabs(CountFunction(func, b));


  if (aX < bX)
  {
    return aX;
  }
  else
  {
    return bX;
  }
}

double CountIterationAlphaFunc(func_t *der, double a, double b)
{
  return 2.0 / FindMaxAbsFunctionBorders(der, a, b);
}

double CountEquivalenceFuncE(func_t *der, double e, double a, double b)
{
  double q = 1.0 - FindMinAbsFunctionBorders(der, a, b) / FindMaxAbsFunctionBorders(der, a, b);
  return e * (1 - q) / q;
}

static double f(double x)
{
  return x;
}

void CreateFunctionIterationFunc(func_t *poly, func_t *iter, double alpha)
{
  CopyFunction(poly, iter);
  MultiplyFunction(iter, -alpha);
  AddElement(iter, &f, 1.0);
}

double FixedPointIteratrionFunc(func_t *func, double e1, double x1, int *iterations)
{
  double x2 = CountFunction(func, x1);
  (*iterations)++;

#ifdef FIRST_ANSWERS
  if (*iterations < 5)
  {
    printf("  x%i = %lf\n", *iterations, x1);
  }
#endif // FIRST_ANSWERS

  if (fabs(x2 - x1) < e1)
  {
    return x2;
  }
  else
  {
    return FixedPointIteratrionFunc(func, e1, x2, iterations);
  }
}

double FixedPointIterationMethod_function(func_t *func, func_t *der, double a, double b, double *e, int *iterations)
{
#ifdef CHECK_INTERVAL
  if (_checkIntervalFunc(func, a, b) == 0 || _checkIntervalDerivativeFunc(der, a, b) == 0)
  {
    *iterations = -1;
    return 0.0;
  }
#endif // CHECK_INTERVAL
  {
    func_t iter;
    double answer;
    double alpha;
    double e1;

    InitFunction(&iter);

    alpha = CountIterationAlphaFunc(der, a, b);
    e1 = CountEquivalenceFuncE(der, *e, a, b);

    if (CountFunction(der, b) < 0.0)
    {
      CreateFunctionIterationFunc(func, &iter, -alpha);
    }
    else
    {
      CreateFunctionIterationFunc(func, &iter, alpha);
    }

    answer = FixedPointIteratrionFunc(&iter, e1, b, iterations);

    DeleteFunction(&iter);

    return answer;
  }
}