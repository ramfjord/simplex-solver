#ifndef STDS 
#define STDS 1
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#endif
#ifndef MATRIX
#define MATRIX 1
#include "matrix.h"
#endif

void print_nxm_to(double *mat, int n, int m, FILE *out){
  for(int r = 0; r < n; r++){
    fputc('[', out);
    for(int c = 0; c < m; c++){
      fprintf(out, "%f",*(mat + r*m + c));
      if(c < m-1)
	fprintf(out, ", ");
    }
    fprintf(out, "]\n");
  }
}

void print_nxm(double *mat, int n, int m){
  print_nxm_to(mat, n, m, stdout);
}

int matrix_equals(double *m1, double *m2, int n, int m){
  for(int i = 0; i < n*m; i++)
    if(m1[i] != m2[i])
      return 0;
  return 1;
}

void add_nxm_matrices(double *m1, double *m2, double *m3, int n, int m){
  for(int r = 0; r < n; r++)
    for(int c = 0; c < m; c++)
      *(m3 + r*m + c) = *(m1 + r*m + c) + *(m2 + r*m + c);
}

void sub_nxm_matrices(double *m1, double *m2, double *m3, int n, int m){
  for(int r = 0; r < n; r++)
    for(int c = 0; c < m; c++)
      *(m3 + r*m + c) = *(m1 + r*m + c) - *(m2 + r*m + c);
}

// sets m1 to the nxm 0 matrix
void zero_matrix(double *m1, int n, int m){
  for(int r = 0; r < n; r++)
    for(int c = 0; c < m; c++)
      *(m1 + r*m + c) = 0;
}

// sets m1 to the nxn identity matrix
void identity_matrix(double *m1, int n){
  zero_matrix(m1,n,n);
  for(int r = 0; r < n; r++)
    *(m1 + r*n + r) = 1;
}

void matrix_copy(double *src, double *dest, int n, int m){
  for(int r = 0; r < n; r++)
    for(int c = 0; c < m; c++)
      *(dest + r*n + c) = *(src + r*n + c);

}

// mult_nxm_by_scalar multiplies a matrix by a scalar
//
// m1 is an n x m matrix
// scalar is the scalar
//
// final values are stored in m1
void mult_nxm_by_scalar(double *m1, double scalar, int n, int m){
  for(int r = 0; r < n; r++)
    for(int c = 0; c < m; c++)
      *(m1 + r*m + c) *= scalar;
}

double dot_product(double *m1, double *m2, int size){
  double ret = 0;
  for(int i = 0; i < size; i++)
    ret += m1[i] * m2[i];
#ifdef PRINT_OPERATIONS/*{{{*/
  printf("in dot_product: %d\n", size);
  print_nxm(m1, 1, size);
  printf("dot\n");
  print_nxm(m2, 1, size);
  printf("= %f\n\n",ret);
#endif/*}}}*/
  return ret;
}

// n is set to the cross product of a and b. a,b,n are all 3 dimensional vectors.
double cross_product_3d(double *a, double *b, double *n){
  n[X_AXIS] = a[Y_AXIS] * b[Z_AXIS] - a[Z_AXIS] * b[Y_AXIS];
  n[Y_AXIS] = a[Z_AXIS] * b[X_AXIS] - a[X_AXIS] * b[Z_AXIS];
  n[Z_AXIS] = a[X_AXIS] * b[Y_AXIS] - a[Y_AXIS] * b[X_AXIS];
}

/*
  scalar = scalar value by which m1 is multiplied
  m1 = 1xn matrix to be multiplied by scalar
  res = scalar*m1 + res
  size = n
*/
void mult_by_scaler_and_add(double scalar, double *m1, double *res, int size){
  for(int i = 0; i < size; i++)
    res[i] += scalar * m1[i];
}

/*
   m1 = 1*n matrix
   m2 = n*n matrix
   m3 = 1*n matrix, where result will be stored
   n = n
*/
void mult_1xn_by_nxn_store(double *m1, double *m2, double *m3, int n){
  for(int c = 0; c < n; c++)
    m3[c] = 0;
  for(int r = 0; r < n; r++)
    mult_by_scaler_and_add(m1[r], (m2 + n*r), m3, n);
#ifdef PRINT_OPERATIONS/*{{{*/
  printf("in 1x%d\n", n);
  print_nxm(m1, 1, n);
  printf("x\n");
  print_nxm(m2, n, n);
  printf("---------\n");
  print_nxm(m3, 1, n);
  putchar('\n');
#endif/*}}}*/
  return;
}

/*
   m1 = n*n matrix
   m2 = n*n matrix
   m3 = n*n matrix, where result will be stored
   n = n
*/
void mult_nxn_store(double *m1, double *m2, double *m3, int n){
#ifdef PRINT_OPERATIONS/*{{{*/
  printf("in %dx%d\n", n, n);
  print_nxm(m1, n, n);
  printf("x\n");
  print_nxm(m2, n, n);
  printf("---------\n");
#endif/*}}}*/
  for(int r = 0; r < n; r++){
    mult_1xn_by_nxn_store((m1 + (n*r)), m2, (m3 + (n*r)), n);
  }
#ifdef PRINT_OPERATIONS/*{{{*/
  print_nxm(m3, n, n);
  putchar('\n');
#endif/*}}}*/
  return;
}

void mult_nxn(double *m1, double *m2, int n){
  double m3[n][n];
  mult_nxn_store(m1, m2, CM m3, n);
  matrix_copy(CM m3, m1, n,n);
}


void mult_2x2_matrix(double m1[2][2], double m2[2][2], double m3[2][2]){
  m3[0][0] = m1[0][0] * m2[0][0] + m1[0][1] * m2[1][0];/*{{{*/
  m3[0][1] = m1[0][0] * m2[0][1] + m1[0][1] * m2[1][1];
  m3[1][0] = m1[1][0] * m2[0][0] + m1[1][1] * m2[1][0];
  m3[1][1] = m1[1][0] * m2[0][1] + m1[1][1] * m2[1][1];
  return;
}/*}}}*/

/*
  Takes a n sized vector *translation, and an (n+1)*(n+1) T where the resulting 
  translation matrix will be stored.
*/
void get_n_translation(double *T, int n, double *translation){ 
  identity_matrix(T,4);
  for(int c = 0; c < n; c++){
    *(T + (n+1)*n + c)= translation[c];
  }
}/*}}}*/

void get_translation(double *T, double x, double y){
  double translation[3];/*{{{*/
  translation[0] = x;
  translation[1] = y;
  translation[2] = 1;
  get_n_translation(T, 2, translation);
}/*}}}*/

// sets T to the (n+1)x(n+1) scaling transition matrix, where each dimension is scaled
// by its corresponding value in scale (an n dimensional vector).
void get_n_scale(double *T, int n, double *scale){ 
  identity_matrix(T,n+1);
  for(int i = 0; i < n; i++){
    *(T + (n+1)*i + i)= scale[i];
  }
}/*}}}*/

void get_scale(double *T, double x, double y){
  double m[2];/*{{{*/
  m[X_AXIS] = x;
  m[Y_AXIS] = y;
  get_n_scale(T, 2, m);
}/*}}}*/

void homogenize(double *V){
  double scale = 1 / V[3];
  mult_nxm_by_scalar(V,scale,1,4);
}
