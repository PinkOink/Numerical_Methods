#ifndef FIXEDPOINTITERATION_FUNC_H_INCLUDED__
#define FIXEDPOINTITERATION_FUNC_H_INCLUDED__
#pragma once
#include "transcendental.h"

double FindMaxAbsFunctionBorders(func_t *func, double a, double b);

double FindMinAbsFunctionBorders(func_t *func, double a, double b);

double CountIterationAlphaFunc(func_t *der, double a, double b);

double CountEquivalenceFuncE(func_t *der, double e, double a, double b);

void CreateFunctionIterationFunc(func_t *poly, func_t *iter, double alpha);

double FixedPointIteratrionFunc(func_t *func, double e1, double x1, int *iterations);

double FixedPointIterationMethod_function(func_t *func, func_t *der, double a, double b, double *e, int *iterations);

#endif