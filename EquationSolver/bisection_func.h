#ifndef BISECTION_FUNC_H_INCLUDED__
#define BISECTION_FUNC_H_INCLUDED__
#pragma once
#include "transcendental.h"

double BisectionMethod_function(func_t *func, double a, double b, double *e, int *iterations);

#endif