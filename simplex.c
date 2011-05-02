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

#include <readline/readline.h>
#include <readline/history.h>

// print table flags
#define INCL_MULT 1
#define NO_Z
#define PHASE_1 2
#define NO_X 4
#define NO_S 8
#define NO_RHS 16
#define TOP_ONLY 32
#define NO_Y 64
#define FIRST_PHASE1 128
#define PIVOT 256

void put_chars(int c, int num){
  for(int i = 0; i < num; i++)
    putchar(c);
}

void print_elt_pretty(double **matrix, int r, int c, int flags){
  if(matrix[r][c] == 0)
    matrix[r][c] = 0;

  if(flags & PIVOT){
    if(matrix[r][c] < 0)
      printf("[%.3f]",matrix[r][c]);
    else
      printf(" [%.3f]",matrix[r][c]);
  }
  else if(matrix[r][c] < 0)
    printf(" %.3f ",matrix[r][c]);
  else
    printf("  %.3f ",matrix[r][c]);
}

void print_mult(double elt){
  if(elt == 0)
    elt = 0;
  if(elt < 0)
    printf("%.3f : ", elt);
  else
    printf(" %.3f : ", elt);
}

void loop(double **matrix, int num_x, int num_s, int phase);

//
// flags should be set to 0. bitwise or it with the options you want.
//    options include INCL_MULT, NO_Z(not implemented), PHASE_1
//
void print_tabl(double **matrix, int num_x, int num_s, int flags, int piv_row, int piv_col){
  int num_cons = num_s + 1;
  int num_tot = num_x + num_s + 3;

  int incl_mult = flags & INCL_MULT;
  int phase_1 = flags & PHASE_1;
  int no_x = flags & NO_X;
  int no_y = flags & NO_Y;
  int no_s = flags & NO_S;
  int top_only = flags & TOP_ONLY;

  if(phase_1) num_x -= 1;

  putchar('\n');

  //
  // print header
  //
  if(incl_mult)
    printf("mult\t");

  printf("    z");

  if(!no_x)
  for(int i = 1; i <= num_x; i++)
    printf("\t |  x%d",i);

  if(!no_y)
  for(int i = 1; i <= num_s; i++)
    printf("\t |  s%d",i);
  
  if(phase_1) printf("\t\ts%d",2); // print x0 if it's phase 1
  puts("\t|   RHS");

  if(top_only) return;

  put_chars('-',8*num_tot);
  putchar('\n');
  
  //
  // print matrix
  //

  for(int r = 1; r < num_cons; r++){
    if(incl_mult)
      print_mult(matrix[r][0]);


    print_elt_pretty(matrix,r,1,0);
    
    if(! no_x);
    for(int c = 2; c < num_x + 2; c++){
      if(r == piv_row && c == piv_col && INCL_MULT)
	print_elt_pretty(matrix, r, c, PIVOT);
      else
	print_elt_pretty(matrix, r, c, 0);
    }
    if(! no_y)
    for(int c = num_x + 2; c < num_tot-1; c++){
      if(r == piv_row && c == piv_col && INCL_MULT)
	print_elt_pretty(matrix, r, c, PIVOT);
      else
	print_elt_pretty(matrix, r, c, 0);
    }

    print_elt_pretty(matrix, r, num_tot-1, 0); // rhs
    putchar('\n');
  }

  put_chars('-',8*num_tot);
  putchar('\n');

  if(incl_mult)
    print_mult(matrix[0][0]);
  for(int c = 1; c < num_tot; c++)
    print_elt_pretty(matrix,0,c,0);

  putchar('\n');
  putchar('\n');
}

void print_t(double **matrix, int num_x, int num_s, int flags){
  print_tabl(matrix, num_x, num_s, flags, (0-1), (0-1));
}

void free_mat(double **mat){
  double *temp;
  while(*mat != NULL){
    temp = *mat;
    mat += 1;
    free(temp);
  }
}
    
int origin_infeasible(double **matrix, int num_x, int num_s){
  int RHS = num_x + num_s + 2;
  for(int i = 1; i < num_s + 1; i++)
    if(matrix[i][RHS] < 0)
      return 1;
  return 0;
}
	  
void next_pivot(double **matrix, int num_x, int num_s, int entering_var, int flags){
  int fphase1 = flags & FIRST_PHASE1;
  int phase_1 = flags & PHASE_1;
  int num_cons = num_s + 1;
  int num_tot = num_x + num_s + 3;
  int piv_col, piv_row;
  double max = 0;

  putchar('\n');

  if(fphase1){ // if this is the first phase of a 2 phase thing
    for(int i = 1; i < num_cons; i++){
      if(matrix[i][num_tot-1] < max){
	max = matrix[i][num_tot-1];
	piv_row = i;
      }
    }
    if(max == 0){
      printf("something smells fishy...\n");
      return;
    }

    piv_col = 2;  // always the row with x_0
  } 

    //
    // Else, we are doing a normal
    //
  else {
  if(origin_infeasible(matrix, num_x, num_s)){
    printf("origin infeasible, starting two-phase solution\n");
    double *matrix_2[100];

    // set two-phase objective function
    matrix_2[0] = (double *)malloc(sizeof(double) * (num_tot+1));
    zero_matrix(matrix_2[0],1,num_tot+1);
    matrix_2[0][1] = -1;
    matrix_2[0][2] = -1;

    //copy the rest of the matrix from the original, adding in x_0
    for(int i = 1; i < num_cons; i++){
      matrix_2[i] = malloc(sizeof(double) * (num_tot+1));
      matrix_copy(matrix[i], matrix_2[i], 1, 2);  // copy mult and z
      matrix_2[i][2] = -1; // add x_0
      matrix_copy(matrix[i]+2, matrix_2[i] + 3, 1, num_tot-2);
    }

    for(int i = num_cons; i < 100; i++)
      matrix_2[i] = NULL;	  // for freeing

    char *buffer;
    next_pivot(matrix_2, num_x+1, num_s, -1, PHASE_1 | FIRST_PHASE1);

    loop(matrix_2, num_x+1, num_s, PHASE_1);

    puts("\n\n###########################");
    puts("completed phase 1 !");
    puts("###########################\n\n");

    //copy matrix 2 back to the original, for phase 2
    for(int i = 1; i < num_cons; i++){
      matrix_copy(matrix_2[i], matrix[i], 1, 2);  // copy mult and z
      matrix_copy(matrix_2[i]+3, matrix[i] + 2, 1, num_tot-2);
    }

    free_mat(matrix_2);
    return;
  }

  if(entering_var == -1){
    for(int i = 2; i < num_tot-1; i++){
      if(matrix[0][i] > max){
	max = matrix[0][i];
	piv_col = i;
      }
    }
    if(max <= 0.0000001){
      puts("but you're done!!!");
      return;
    }
  } else piv_col = entering_var;

  max = 0;
  for(int i = 1; i < num_cons; i++){
     double ratio = matrix[i][piv_col] / matrix[i][num_tot - 1];
     if(ratio > max){
       piv_row = i;
       max = ratio;
     }
  }
  if(max <= 0.00000001){
    puts("no leaving variable! solution is UNBOUNDED!");
    int vchar = (piv_col >= num_x + 2) ? 's' : 'x';
    int subsc = (piv_col >= num_x + 2) ? piv_col-num_x-1 : piv_col -1;
    printf("we can increase %c%d indefinitely\n", vchar, subsc);
    return;
  }

  } // if origin is feasible
    
  printf("pivoting on row %d, col %d, ratio = %.3f / %.3f = %.3f : %.3f\n", piv_row, piv_col, matrix[piv_row][piv_col], matrix[piv_row][num_tot-1],matrix[piv_row][piv_col]/matrix[piv_row][num_tot-1], max);

  //
  // get the multiplier on each col
  //
  for(int i = 0; i < num_cons; i++){
    if(i == piv_row)
      matrix[piv_row][0] = 1.00 / matrix[piv_row][piv_col];
    else
      matrix[i][0] = -matrix[i][piv_col] / matrix[piv_row][piv_col];
  }
  
  print_tabl(matrix, num_x, num_s, INCL_MULT, piv_row, piv_col);
  

  //
  // compute the new matrix
  //
  for(int r = 0; r < num_cons; r++){
    if(r == piv_row)
      continue;
    else
      for(int c = 1; c < num_tot; c++){
	matrix[r][c] += matrix[r][0]*matrix[piv_row][c];
      }
    matrix[r][0] = 0;
  }

  for(int c = 1; c < num_tot; c++){
    matrix[piv_row][c] *= matrix[piv_row][0];
  }
  matrix[piv_row][0] = 0;

  
}

void loop(double **matrix, int num_x, int num_s, int phase){
  if(phase == 1 || phase == PHASE_1)
    phase = PHASE_1;
  else
    phase = 0;
  char *buffer;
  while(1){
    buffer = (phase == PHASE_1) ? readline("compute phase 1 $ ") : readline("compute $ ");
    add_history(buffer);
    if(strstr(buffer, "quit"))
      break;
    if(strstr(buffer, "next"))
      next_pivot(matrix, num_x, num_s, -1, phase);
    if(strstr(buffer, "ent")){
      double entering_var;
      if(1 == sscanf(buffer, "ent %lf", &entering_var))
	next_pivot(matrix, num_x, num_s, entering_var, phase);
    }

    if(strstr(buffer, "help")){
      printf("commands are: \n");
      printf(" quit       : quit's the current phase\n");
      printf(" next       : executes a pivot on the logical entering variable\n");
      printf(" ent <num>  : executes a pivot on column <num>.  note that column 1 is z, and the last column is RHS\n");
      printf(" print      : prints out the current matrix\n");
    }

    if(strstr(buffer, "print")){
      print_t(matrix, num_x, num_s, 0);
    }
    free(buffer);
  }
}

int main(){
  char *buffer;
  while(1){
    if(!(buffer = readline("Simplex$ ")))
      break;
    add_history(buffer);

    if(strstr(buffer, "help")){
      printf("possible commands are:\n");
      puts(" quit          : exit the program");
      puts(" start <v> <c> : start simplex method with <v> variables and <c> constraints");
    }

    if(strstr(buffer, "quit"))
      exit(0);
    
    if(strstr(buffer, "start ")){
      int num_x, num_s, num_tot;
      int num_cons = 0;
      double *matrix[100];
      for(int i = 0; i < 100; i++) matrix[i] = (double *) NULL;

      if(2 != sscanf(buffer, "start %d %d", &num_x, &num_s)){
	puts("usage: $ start <num_x> <num_s>");
	continue;
      }
      num_tot = num_x + num_s + 3; // total number of variables; 3 for multiplier, z, and RHS

      //
      // Read in objective function
      //
      putchar('\n');
      printf("input objective function: z = ");
      for(int i = 0; i < num_x; i++)
	printf(" a%d*x%d",i,i);
      putchar('\n');

      puts("enter the following separated by spaces:");
      printf("           ");
      for(int i = 0; i < num_x; i++)
	printf(" a%d",i);
      puts("");

      buffer = readline("input_args$ ");
      add_history(buffer);

      char *temp = buffer;
      matrix[num_cons] = (double *)malloc(sizeof(double) * num_tot);
      zero_matrix(matrix[num_cons], 1, num_tot);
      for(int i = 1; i <= num_x; i++)
	matrix[num_cons][i+1] = strtod(temp, &temp);
      matrix[num_cons][1] = -1;
      num_cons++;

      //
      // Read in constraints
      //
      putchar('\n');
      printf("the constraints are in the form: s_i");
      for(int i = 0; i < num_x; i++)
	printf(" + a_%dx%d",i,i);
      puts(" = RHS");
      puts("enter the following separated by spaces:");
      printf("     ");
      for(int i = 0; i < num_x; i++)
	printf(" a_%d",i,i);
      puts(" RHS");
      free(buffer);

      for(int which_con = 0; which_con < num_s; which_con++){
	char prompt[10];
	sprintf(prompt,"s%d : ",num_cons);
	buffer = readline(prompt);
	add_history(buffer);

	matrix[num_cons] = (double *)malloc(sizeof(double) * num_tot);
	zero_matrix(matrix[num_cons], 1, num_tot);
	temp = buffer;
	for(int i = 2; i < num_x+2; i++)
	  matrix[num_cons][i] = strtod(temp, &temp);
	matrix[num_cons][num_tot-1] = strtod(temp, &temp); // set rhs
	matrix[num_cons][num_x + 2 + which_con] = 1; // set slack variable to 1

	num_cons++;
	free(buffer);
      }
      for(int i = num_cons; i < 100; i++)
	matrix[i] = NULL;	  // for freeing

      puts("ready to compute: type next to iterate, quit to quit, and finish to finish");

      loop(matrix, num_x, num_s, 0);

      free_mat(matrix);
    } // if start

  } // while
} // main
