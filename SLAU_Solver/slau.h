#ifndef SLAU_H_INCLUDED__
#define SLAU_H_INCLUDED__
#pragma once

#include "matrix_vector.h"
#include <stdio.h>

//matrix_t InitTriagMatrix(void);

matrix_t GetMatrix(FILE *f);

vector_t InitRandVector(void);

void PrepareMatrix(matrix_t *triag);

matrix_t InitRotationMatrix(int i, int j, double tg);

vector_t DeltaVector(vector_t vec);

matrix_t DeltaMatrix(matrix_t mat);

matrix_t CreateRotationMatrix(matrix_t mat);

vector_t SolveSLAU(matrix_t mat, matrix_t rotat, vector_t b);

#endif