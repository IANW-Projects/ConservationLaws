// This project is licensed under the terms of the Creative Commons CC BY-NC-ND 4.0 license.

/* Contains the coefficients of the fourth order diagonal norm derivative operator of
@article{mattsson2004summation,
  title={Summation by parts operators for finite difference approximations of
         second derivatives},
  author={Mattsson, Ken and Nordstr{\"o}m, Jan},
  journal={Journal of Computational Physics},
  volume={199},
  number={2},
  pages={503--540},
  year={2004},
  publisher={Elsevier},
  doi={10.1016/j.jcp.2004.03.001}
}
*/

//------------------------------------------------------------------------------
// Periodic Derivative Operator
//------------------------------------------------------------------------------
#define ORDER 4
#define NUM_BOUNDS_P 0
#define STENCIL_WIDTH_P 5

/*
Coefficients of the first derivative operator.

For each row, the central coefficient corresponds to the nodes at which the
derivative is computed.
*/
REAL constant SBP_diff_P[2*NUM_BOUNDS_P+1][STENCIL_WIDTH_P] = {
  // central coefficients
  {1.0/12.0, -2.0/3.0, 0.0, 2.0/3.0, -1.0/12.0},
};

/*
Coefficients of the inverse mass/norm matrix.


The central coefficient corresponds to the inner nodes of the grid and is 1.
*/
REAL constant M_INV_P[2*NUM_BOUNDS_P+1] = {
  // central coefficients
  1.0
};


//------------------------------------------------------------------------------
// Non-Periodic Derivative Operator
//------------------------------------------------------------------------------

#define ORDER 4
#define NUM_BOUNDS_B 4
#define STENCIL_WIDTH_B 7

/*
Coefficients of the first derivative operator.

For each row, the central coefficient corresponds to the nodes at which the
derivative is computed.
*/
REAL constant SBP_diff_B[2*NUM_BOUNDS_B+1][STENCIL_WIDTH_B] = {
  // left boundary coefficients
  {0.0, 0.0, 0.0, -24.0/17.0, 59.0/34.0, -4.0/17.0, -3.0/34.0},
  {0.0, 0.0, -1.0/2.0, 0.0, 1.0/2.0, 0.0, 0.0},
  {0.0, 4.0/43.0, -59.0/86.0, 0.0, 59.0/86.0, -4.0/43.0, 0.0},
  {3.0/98.0, 0.0, -59.0/98.0, 0.0, 32.0/49.0, -4.0/49.0, 0.0},
  // central coefficients
  {0.0, 1.0/12.0, -2.0/3.0, 0.0, 2.0/3.0, -1.0/12.0, 0.0},
  // right boundary coefficients
  {0.0, 4.0/49.0, -32.0/49.0, 0.0, 59.0/98.0, 0.0, -3.0/98.0},
  {0.0, 4.0/43.0, -59.0/86.0, 0.0, 59.0/86.0, -4.0/43.0, 0.0},
  {0.0, 0.0, -1.0/2.0, 0.0, 1.0/2.0, 0.0, 0.0},
  {3.0/34.0, 4.0/17.0, -59.0/34.0, 24.0/17.0, 0.0, 0.0, 0.0}
};

/*
Coefficients of the inverse mass/norm matrix.

The central coefficient corresponds to the inner nodes of the grid and is 1.
*/
REAL constant M_INV_B[2*NUM_BOUNDS_B+1] = {
  // left boundary coefficients
  48.0/17.0,
  48.0/59.0,
  48.0/43.0,
  48.0/49.0,
  // central coefficients
  1.0,
  // right boundary coefficients
  48.0/49.0,
  48.0/43.0,
  48.0/59.0,
  48.0/17.0
};


//------------------------------------------------------------------------------
// Coeffcients of the high-order dissipation operator
//------------------------------------------------------------------------------

#define STENCIL_WIDTH_HOD 3
#define NUM_BOUNDS_HOD 4


REAL constant D_HO[2*NUM_BOUNDS_HOD+1][STENCIL_WIDTH_HOD] = {
  // left boundary coefficients
  { 0.0, 0.0, -1.0},
  { 0.0, 2.0, -1.0},
  {-1.0, 2.0, -1.0},
  {-1.0, 2.0, -1.0},
  // central coefficients
  {-1.0, 2.0, -1.0},
  // right boundary coefficients
  {-1.0, 2.0, -1.0},
  {-1.0, 2.0, -1.0},
  {-1.0, 2.0,  0.0},
  {-1.0, 0.0,  0.0}
};

