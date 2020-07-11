#include "matrix_vector.h"
#include "slau.h"
#include <stdio.h>
#pragma warning(disable:4996)
//#define TESTING

void MakeMatrix(matrix_t mat)
{
  mat[0][0] = 1;
  mat[0][1] = 2;
  mat[0][2] = -1;
  mat[1][0] = 2;
  mat[1][1] = -3;
  mat[1][2] = 2;
  mat[2][0] = 3;
  mat[2][1] = 1;
  mat[2][2] = 1;
}

void MakeVector(vector_t vec)
{
  vec[0] = 2;
  vec[1] = 2;
  vec[2] = 8;
}

void Solving(void)
{
  matrix_t A = InitIdentityMatrix();
  matrix_t Q;
  vector_t B = InitVector();
  vector_t X;
  vector_t D;
  double normX, normA, normB;

#ifdef TESTING
  MakeMatrix(A);
  MakeVector(B);

#else
  FILE *f = fopen("10 10000000000.txt", "r");

  DeleteMatrix(A);
  DeleteVector(B);
  A = GetMatrix(f);
  B = InitRandVector();
  PrepareMatrix(&A);
//  PrintMatrix(A);
//  PrintVector(B);
  fclose(f);
#endif
//  PrintMatrix(A);
  Q = CreateRotationMatrix(A);
  X = SolveSLAU(A, Q, B);
  D = Matrix_x_Vector(A, X);
  Vector_x_Number(B, -1);
  Vector_plus_Vector(D, B);
  Vector_x_Number(B, -1);

  normA = InfiniteNormMatrix(A);
  normB = InfiniteNormVector(B);
  normX = InfiniteNormVector(X);

  printf("\tI A * X = B; Size = %i\n\n", SIZE);
  printf("D = A * X - B\nNorm(D) = %lf\n\n", InfiniteNormVector(D));
  printf("||A|| = %lf\n", normA);
  printf("||B|| = %lf\n", normB);
  printf("||X|| = %lf\n\n", normX);

  {
    vector_t B1 = InitVector();
    vector_t deltaB1;
    vector_t X1;
    vector_t deltaX1 = InitVector();
    double normDeltaX1, normDeltaB1;

    CopyVectors(B, B1);
    deltaB1 = DeltaVector(B1);
    Vector_plus_Vector(B1, deltaB1);

    X1 = SolveSLAU(A, Q, B1);
    CopyVectors(X1, deltaX1);
    Vector_x_Number(X, -1);
    Vector_plus_Vector(deltaX1, X);
    Vector_x_Number(X, -1);

    normDeltaX1 = InfiniteNormVector(deltaX1);
    normDeltaB1 = InfiniteNormVector(deltaB1);

    printf("\tII A * X = (B + deltaB)\n\n");
    printf("||X|| = %lf\n", normX);
    printf("||deltaX|| = %lf\n", normDeltaX1);
    printf("||B|| = %lf\n", normB);
    printf("||deltaB|| = %lf\n", normDeltaB1);
    printf("k1 = %lf\n\n", normDeltaX1 * normB / (normX * normDeltaB1));

    DeleteVector(X1);
    DeleteVector(deltaX1);
    DeleteVector(B1);
    DeleteVector(deltaB1);
  }

  {
    matrix_t A2 = InitMatrix();
    matrix_t deltaA2;
    matrix_t Q2;
    vector_t X2;
    vector_t deltaX2 = InitVector();
    double normA2, normDeltaA2, normDeltaX2;

    CopyMatrix(A, A2);

    deltaA2 = DeltaMatrix(A2);
    Matrix_plus_Matrix(A2, deltaA2);

    Q2 = CreateRotationMatrix(A2);
    X2 = SolveSLAU(A2, Q2, B);
    CopyVectors(X2, deltaX2);
    Vector_x_Number(X, -1);
    Vector_plus_Vector(deltaX2, X);
    Vector_x_Number(X, -1);

    normA2 = InfiniteNormMatrix(A2);
    normDeltaA2 = InfiniteNormMatrix(deltaA2);
    normDeltaX2 = InfiniteNormVector(deltaX2);

    printf("\tIII (A + deltaA) * X = B\n\n");
    printf("||X|| = %lf\n", normX);
    printf("||deltaX|| = %lf\n", normDeltaX2);
    printf("||deltaA|| = %lf\n", normDeltaA2);
    printf("||deltaA + A|| = %lf\n", normA2);
    printf("k2 = %lf\n\n", normDeltaX2 * normA2 / (normX * normDeltaA2));

    DeleteMatrix(A2);
    DeleteMatrix(deltaA2);
    DeleteMatrix(Q2);
    DeleteVector(X2);
    DeleteVector(deltaX2);
  }

  DeleteMatrix(A);
  DeleteMatrix(Q);
  DeleteVector(B);
  DeleteVector(X);
}

int main(void)
{
  Solving();
  return 0;
}