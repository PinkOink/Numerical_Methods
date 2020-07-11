#include "polynomial.h"
#include "bisection_poly.h"
#include <math.h>
#include <stdio.h>

static int _checkNegativeCoefExist(polynomial_t *polynomial)
{
  monomial_t *buf = polynomial->first;

  while (buf != polynomial->last)
  {
    if (buf->coef < 0.0)
    {
      return 1;
    }
    buf = buf->next;
  }

  return 0;
}

static double _searchMaxNegativeCoef(polynomial_t *polynomial)
{
  monomial_t *buf = polynomial->first;
  double max = 0.0;

  while (buf != polynomial->last)
  {
    if (buf->coef < max)
    {
      max = buf->coef;
    }
    buf = buf->next;
  }

  return -max;
}

static int _searchFirstNegativeCoef(polynomial_t *polynomial)
{
  monomial_t *buf = polynomial->first;
  int i = 0;

  while (buf != polynomial->last)
  {
    if (buf->coef < 0.0)
    {
      return i;
    }
    i++;
    buf = buf->next;
  }

  return -1;
}

static void _changeSignX(polynomial_t *polynomial)
{
  monomial_t *buf = polynomial->first;

  while (buf != polynomial->last)
  {
    if (buf->pow % 2 == 1)
    {
      buf->coef = -buf->coef;
    }
    buf = buf->next;
  }
}

double SearchUpperBoundPositive(polynomial_t *polynomial)
{
  if (polynomial->first->coef < 0.0)
  {
    MultiplyPolinomial(polynomial, -1.0);
  }

  if (_checkNegativeCoefExist(polynomial))
  {
    int m = _searchFirstNegativeCoef(polynomial);
    double a0 = polynomial->first->coef;
    double a = _searchMaxNegativeCoef(polynomial);

    return (1 + pow(a / a0, 1.0 / m));
  }
  else
  {
    return 0.0;
  }
}

double SearchLowerBoundPositive(polynomial_t *polynomial)
{
  polynomial_t reverse;
  double bound;

  InitPolynomial(&reverse);
  ReversePolynomial(polynomial, &reverse);
  bound = SearchUpperBoundPositive(&reverse);
  DeletePolynomial(&reverse);

  if (bound != 0.0)
  {
    return 1 / bound;
  }
  else
  {
    return 0.0;
  }
}

void FindPositiveBounds(polynomial_t *poly, double *a, double *b)
{
  *a = SearchLowerBoundPositive(poly);
  *b = SearchUpperBoundPositive(poly);
}

double SearchUpperBoundNegative(polynomial_t *polynomial)
{
  double bound;

  _changeSignX(polynomial);
  bound = SearchLowerBoundPositive(polynomial);
  _changeSignX(polynomial);

  return -bound;
}

double SearchLowerBoundNegative(polynomial_t *polynomial)
{
  double bound;

  _changeSignX(polynomial);
  bound = SearchUpperBoundPositive(polynomial);
  _changeSignX(polynomial);

  return -bound;
}

void FindNegativeBounds(polynomial_t *poly, double *a, double *b)
{
  *a = SearchLowerBoundNegative(poly);
  *b = SearchUpperBoundNegative(poly);
}

double BisectionMethod_polynomial(polynomial_t *polynomial, double a, double b, double *e, int *iterations)
{
  double c = (b + a) / 2;

  (*iterations)++;

  if (fabs(b - a) < 2 * *e || CountPolynomial(polynomial, c) == 0)
  {
    return c;
  }
  else
  {
    if (CountPolynomial(polynomial, c) * CountPolynomial(polynomial, a) < 0.0)
    {
      return BisectionMethod_polynomial(polynomial, a, c, e, iterations);
    }
    else
    {
      return BisectionMethod_polynomial(polynomial, c, b, e, iterations);
    }
  }
}