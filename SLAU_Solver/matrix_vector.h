#ifndef MATRIX_VECTOR_H_INCLUDED__
#define MATRIX_VECTOR_H_INCLUDED__
#pragma once

#define SIZE 10

typedef double *vector_t;

typedef double **matrix_t;

vector_t InitVector(void);

matrix_t InitMatrix(void);

matrix_t InitIdentityMatrix(void);

void CopyVectors(vector_t from, vector_t to);

void CopyMatrix(matrix_t from, matrix_t to);

void PrintVector(vector_t vec);

void PrintMatrix(matrix_t mat);

void TransposeMatrix(matrix_t mat);

void Vector_plus_Vector(vector_t A, vector_t B);

void Matrix_plus_Matrix(matrix_t A, matrix_t B);

void Vector_x_Number(vector_t vec, double n);

void Matrix_x_Number(matrix_t mat, double n);

matrix_t Vector_x_VectorTRANSPOSE(vector_t vecA, vector_t vecB);

vector_t Matrix_x_Vector(matrix_t mat, vector_t vec);

matrix_t Matrix_x_Matrix(matrix_t matA, matrix_t matB);

double EuclideanNormVector(vector_t vec);

double EuclideanNormMatrix(matrix_t diag);

double InfiniteNormVector(vector_t vec);

double InfiniteNormMatrix(matrix_t mat);

void DeleteVector(vector_t vec);

void DeleteMatrix(matrix_t mat);

#endif