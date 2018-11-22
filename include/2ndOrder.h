// This project is licensed under the terms of the Creative Commons CC BY-NC-ND 4.0 license.

// Contains the coefficients of second order derivative operator

//------------------------------------------------------------------------------
// Derivative Operator
//------------------------------------------------------------------------------

#define ORDER 2
#define NUM_BOUNDS 1
#define STENCIL_WIDTH 3

/*
Coefficients of the first derivative operator.

For each row, the central coefficient corresponds to the nodes at which the
derivative is computed.
*/
REAL constant SBP_diff[2*NUM_BOUNDS+1][STENCIL_WIDTH] = {
  // left boundary coefficients
  {0.0, -1.0, 1.0},
  // central coefficients
  {-1.0/2.0, 0.0, 1.0/2.0},
  // right boundary coefficients
  {-1.0, 1.0, 0.0}
};


/*
Coefficients of the inverse mass/norm matrix.

The central coefficient corresponds to the inner nodes of the grid and is 1.
*/
REAL constant M_INV[2*NUM_BOUNDS+1] = {
  // left boundary coefficients
  2.0,
  // central coefficients
  1.0,
  // right boundary coefficients
  2.0
};
