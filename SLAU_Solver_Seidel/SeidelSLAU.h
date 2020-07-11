#ifndef SEIDELSLAU_H_INCLUDED__
#define SEIDELSLAU_H_INCLUDED__
#pragma once
#include "matrix.h"

vector_t InitRandVector(void);

int CheckDiagDominion(matrix_t mat);

matrix_t CreateIterationMatrixD(matrix_t mat);

vector_t CreateGVectorD(matrix_t mat, vector_t b);

matrix_t CreateIterationMatrix(matrix_t mat, double a);

vector_t CreateGVector(vector_t b, double a);

void SeidelIteration(vector_t x0, matrix_t C, vector_t g, double e0, int *iterations);

void PrepareMatrix(matrix_t *triag);

#endif