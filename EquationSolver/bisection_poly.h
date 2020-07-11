#ifndef BISECTION_POLY_H_INCLUDED__
#define BISECTION_POLY_H_INCLUDED__
#pragma once
#include "polynomial.h"

double SearchUpperBoundPositive(polynomial_t *polynomial);

double SearchLowerBoundPositive(polynomial_t *polynomial);

void FindPositiveBounds(polynomial_t *poly, double *a, double *b);

double SearchUpperBoundNegative(polynomial_t *polynomial);

double SearchLowerBoundNegative(polynomial_t *polynomial);

void FindNegativeBounds(polynomial_t *poly, double *a, double *b);

double BisectionMethod_polynomial(polynomial_t *polynomial, double a, double b, double *e, int *iterations);

#endif