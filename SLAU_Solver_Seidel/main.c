#include "matrix.h"
#include "SeidelSLAU.h"
#define E0(e, norm) (e * (1.0 - norm) / norm)
#define det     0.0010
#pragma warning(disable: 4996)

void main(void)
{
  matrix_t A;
  matrix_t Atrans = InitMatrix();
  vector_t b;
  vector_t x = InitVector();

  matrix_t buf;
  vector_t bbuf = InitRandVector();

  matrix_t L, R, D, DInv;

  FILE *f = fopen("in.txt", "r");

  matrix_t C;
  vector_t g = InitVector();

  int iterations = 0;
  double e = 0.001;
//  double normCinf;

  buf = GetMatrix(f);
  Matrix_x_Number(buf, 0.1);
  CopyMatrix(buf, Atrans);
  TransposeMatrix(Atrans);

  A = Matrix_x_Matrix(Atrans, buf);
  b = Matrix_x_Vector(Atrans, bbuf);
  TransposeMatrix(Atrans);
  DeleteMatrix(buf);

  L = DownTriagMatrix(A);
  R = UpTriagMatrix(A);
  D = DiagMatrix(A);
  DInv = InverseDiagMatrix(D);

  buf = Matrix_plus_Matrix(L, R);
  C = Matrix_x_Matrix(DInv, buf);
  Matrix_x_Number(C, -1);
  g = Matrix_x_Vector(DInv, b);
  CopyVectors(g, x);

  SeidelIteration(x, C, g, e, &iterations);

  {
    vector_t a = Matrix_x_Vector(Atrans, x);
    vector_t P;

    Vector_x_Number(a, -1);
    P = Vector_plus_Vector(a, bbuf);

    printf("det(A) = %g\n", det);
    printf("||P|| = %lf\n", InfiniteNormVector(P));
    printf("Iterations = %i\n", iterations);

    DeleteVector(a);
    DeleteVector(P);
  }

  fclose(f);
  DeleteMatrix(A);
  DeleteMatrix(Atrans);
  DeleteVector(bbuf);
  DeleteVector(b);
  DeleteVector(x);
  DeleteMatrix(C);
  DeleteVector(g);
  DeleteMatrix(L);
  DeleteMatrix(R);
  DeleteMatrix(D);
  DeleteMatrix(DInv);
}