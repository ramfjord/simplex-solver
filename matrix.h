#define CM (double *)

#define X_AXIS 0
#define Y_AXIS 1
#define Z_AXIS 2
//#define PRINT_OPERATIONS

#define W_LEFT 0
#define W_RIGHT 1
#define W_TOP 2
#define W_BOT 3
#define W_FRONT 4
#define W_BACK 5




//
// Headers
//

// print matrices
void print_nxm(double *mat, int n, int m);
void print_nxm_to(double *mat, int n, int m, FILE *out);

//matrix operations

// return 1 if m1 = m2, 0 otherwise
int matrix_equals(double *m1, double *m2, int n, int m);

// stores m1 + m2 in m3, where each is an n*m matrix
void add_nxm_matrices(double *m1, double *m2, double *m3, int n, int m);

// stores m1 - m2 in m3, where each is an n*m matrix
void sub_nxm_matrices(double *m1, double *m2, double *m3, int n, int m);

// multiplies m1, an nxm matrix, by scalar
void mult_nxm_by_scalar(double *m1, double scalar, int n, int m);

//returns the dot product of m1, m2, where each is a 1xsize matrix
double dot_product(double *m1, double *m2, int size);

/*
   multiplies m1 by m2 to get m3
   m1,m3 are 1*n matrices.
   m2 is an n*n matrix
*/
void mult_1xn_by_nxn_store(double *m1, double *m2, double *m3, int n);

/*
   multiplies m1 by m2 to get m3
   m1,m2,m3 are all n*n matrices
*/
void mult_nxn_store(double *m1, double *m2, double *m3, int n);
void mult_2x2_matrix(double m1[2][2], double m2[2][2], double m3[2][2]);

// converts n*m matrix m1 to the n*m zero matrix
void zero_matrix(double *m1, int n, int m);

// converts n*n matrix m1 to the n*n identity matrix
void identity_matrix(double *m1, int n);

// copies src to dest, where src and dest are nxm matrices
void matrix_copy(double *scr, double *dest, int n, int m);

void get_viewer_transition_matrix(double *T, double *obj_loc, double *viewer_orig);

