#include "fixedpointiteration_poly.h"
#include "polynomial.h"
#include "bisection_poly.h"
#include <math.h>
#include <stdio.h>
#define CHECK_INTERVAL
#define FIRST_ANSWERS

static int _checkIntervalPoly(polynomial_t *poly, double a, double b)
{
  if (CountPolynomial(poly, a) * CountPolynomial(poly, b) < 0.0)
  {
    return 1;
  }
  else
  {
    return 0;
  }
}

static int _checkIntervalDerivativePoly(polynomial_t *poly, double a, double b)
{
  polynomial_t der;
  double delta = 0.01;

  InitPolynomial(&der);
  DerivativePolynomial(poly, &der);
  if (CountPolynomial(&der, a) * CountPolynomial(&der, b) < 0.0)
  {
    DeletePolynomial(&der);
    return 0;
  }
  else
  {
    while (a + delta < b)
    {
      if (CountPolynomial(&der, a) * CountPolynomial(&der, a + delta) < 0.0)
      {
        DeletePolynomial(&der);
        return 0;
      }
      a += delta;
    }
    DeletePolynomial(&der);
    return 1;
  }
}

double FindMaxAbsPolynomialBorders(polynomial_t *poly, double a, double b)
{
  double aX = fabs(CountPolynomial(poly, a));
  double bX = fabs(CountPolynomial(poly, b));

  if (aX > bX)
  {
    return aX;
  }
  else
  {
    return bX;
  }
}

double FindMinAbsPolynomialBorders(polynomial_t *poly, double a, double b)
{
  double aX = fabs(CountPolynomial(poly, a));
  double bX = fabs(CountPolynomial(poly, b));


  if (aX < bX)
  {
    return aX;
  }
  else
  {
    return bX;
  }
}

double CountIterationAlphaPoly(polynomial_t *der, double a, double b)
{
  return 2.0 / FindMaxAbsPolynomialBorders(der, a, b);
}

double CountEquivalencePolyE(polynomial_t *der, double e, double a, double b)
{
  double q = 1.0 - FindMinAbsPolynomialBorders(der, a, b) / FindMaxAbsPolynomialBorders(der, a, b);
  return e * (1 - q) / q;
}

void CreateFunctionIterationPoly(polynomial_t *poly, polynomial_t *iter, double alpha)
{
  CopyPolynomial(poly, iter);
  MultiplyPolinomial(iter, -alpha);
  AddMonomial(iter, 1.0, 1);
}

double FixedPointIteratrionPoly(polynomial_t *func, double e1, double x1, int *iterations)
{
  double x2 = CountPolynomial(func, x1);
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
    return FixedPointIteratrionPoly(func, e1, x2, iterations);
  }
}

double FixedPointIterationMethod_polynomial(polynomial_t *poly, double a, double b, double *e, int *iterations)
{
#ifdef CHECK_INTERVAL
  if ((_checkIntervalPoly(poly, a, b) == 0 || _checkIntervalDerivativePoly(poly, a, b) == 0) && *iterations == 0)
  {
    *iterations = -1;
    return 0.0;
  }
#endif // CHECK_INTERVAL
  {
    polynomial_t der;
    polynomial_t iter;
    double answer;
    double alpha;
    double e1;

    InitPolynomial(&der);
    InitPolynomial(&iter);
    DerivativePolynomial(poly, &der);
    printf("%f %f\n", a, b);
    alpha = CountIterationAlphaPoly(&der, a, b);
    e1 = CountEquivalencePolyE(&der, *e, a, b);

    if (CountPolynomial(&der, b) < 0.0)
    {
      CreateFunctionIterationPoly(poly, &iter, -alpha);
    }
    else
    {
      CreateFunctionIterationPoly(poly, &iter, alpha);
    }
    printf("%f\n", e1);
    answer = FixedPointIteratrionPoly(&iter, e1, b, iterations);

    DeletePolynomial(&der);
    DeletePolynomial(&iter);

    return answer;
  }
}