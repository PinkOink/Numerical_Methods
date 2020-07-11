#ifndef FIXEDPOINTITERATION_POLY_H_INCLUDED__
#define FIXEDPOINTITERATION_POLY_H_INCLUDED__
#pragma once
#include "polynomial.h"

double FindMaxAbsPolynomialBorders(polynomial_t *poly, double a, double b);

double FindMinAbsPolynomialBorders(polynomial_t *poly, double a, double b);

double CountIterationAlphaPoly(polynomial_t *der, double a, double b);

double CountEquivalencePolyE(polynomial_t *der, double e, double a, double b);

/*      func(x) = x - a * poly(x)      */
/*      sign (a) = sign (polyDer(x))      */
void CreateFunctionIterationPoly(polynomial_t *poly, polynomial_t *iter, double alpha);

double FixedPointIteratrionPoly(polynomial_t *func, double e1, double x1, int *iterations);

double FixedPointIterationMethod_polynomial(polynomial_t *poly, double a, double b, double *e, int *iterations);

#endif