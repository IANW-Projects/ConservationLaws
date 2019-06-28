// This project is licensed under the terms of the Creative Commons CC BY-NC-ND 4.0 license.

/* Contains the coefficients of the extended second order diagonal norm derivative operator of
@article{mattsson2014optimal,
  title={Optimal diagonal-norm {SBP} operators},
  author={Mattsson, Ken and Almquist, Martin and Carpenter, Mark H},
  journal={Journal of Computational Physics},
  volume={264},
  pages={91--111},
  year={2014},
  publisher={Elsevier},
  doi={10.1016/j.jcp.2013.12.041}
}
*/


//------------------------------------------------------------------------------
// Periodic Derivative Operator
//------------------------------------------------------------------------------
#define ORDER 2
#define NUM_BOUNDS_P 0
#define STENCIL_WIDTH_P 3

/*
Coefficients of the first derivative operator.

For each row, the central coefficient corresponds to the nodes at which the
derivative is computed.
*/
REAL constant SBP_diff_P[2*NUM_BOUNDS_P+1][STENCIL_WIDTH_P] = {
  // central coefficients
  {-1.0/2.0, 0.0, 1.0/2.0},
};


/*
Coefficients of the inverse mass/norm matrix.


The central coefficient corresponds to the inner nodes of the grid and is 1.
*/
REAL constant M_INV_P[2*NUM_BOUNDS_P+1] = {
  // central coefficients
  1.0,
};


//------------------------------------------------------------------------------
// Non-Periodic Derivative Operator
//------------------------------------------------------------------------------

#define ORDER 2
#define NUM_BOUNDS_B 3
#define STENCIL_WIDTH_B 5

/*
Coefficients of the first derivative operator.

For each row, the central coefficient corresponds to the nodes at which the
derivative is computed.
*/
REAL constant SBP_diff_B[2*NUM_BOUNDS_B+1][STENCIL_WIDTH_B] = {
  // left boundary coefficients
  {0.0, 0.0, -6.0/5.0, 7.0/5.0, -1.0/5.0},
  {0.0, -1.0/2.0, 0.0, 1.0/2.0, 0.0},
  {1.0/11.0, -7.0/11.0, 0.0, 6.0/11.0, 0.0},
  // central coefficients
  {0.0, -1.0/2.0, 0.0, 1.0/2.0, 0.0},
  // right boundary coefficients
  {0.0, -6.0/11.0, 0.0, 7.0/11.0, -1.0/11.0},
  {0.0, -1.0/2.0, 0.0, 1.0/2.0, 0.0},
  {1.0/5.0, -7.0/5.0, 6.0/5.0, 0.0, 0.0}
};


/*
Coefficients of the inverse mass/norm matrix.

The central coefficient corresponds to the inner nodes of the grid and is 1.
*/
REAL constant M_INV_B[2*NUM_BOUNDS_B+1] = {
  // left boundary coefficients
  12.0/5.0,
  6.0/7.0,
  12.0/11.0,
  // central coefficients
  1.0,
  // right boundary coefficients
  12.0/11.0,
  6.0/7.0,
  12.0/5.0
};


//------------------------------------------------------------------------------
// Dissipation Operator
//------------------------------------------------------------------------------

#define NUM_BOUNDS_HOD 3
#define STENCIL_WIDTH_HOD 2

REAL constant D_HO[2*NUM_BOUNDS_HOD+1][STENCIL_WIDTH_HOD] = {
  // left boundary coefficients
  { 0.0, 1.0},
  {-1.0, 1.0},
  {-1.0, 1.0},
  // central coefficients
  {-1.0, 1.0},
  // right boundary coefficients
  {-1.0, 1.0},
  {-1.0, 1.0},
  {-1.0, 0.0}
};

