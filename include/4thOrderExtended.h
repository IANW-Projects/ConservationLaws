// This project is licensed under the terms of the Creative Commons CC BY-NC-ND 4.0 license.

/* Contains the coefficients of the extended fourth order diagonal norm derivative operator of
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
  {0.08333333333333333, -0.6666666666666666, 0.0, 0.6666666666666666, -0.08333333333333333},
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

#define ORDER 4
#define NUM_BOUNDS_B 6
#define STENCIL_WIDTH_B 11

/*
Coefficients of the first derivative operator.

For each row, the central coefficient corresponds to the nodes at which the
derivative is computed.
*/
REAL constant SBP_diff_B[2*NUM_BOUNDS_B+1][STENCIL_WIDTH_B] = {
  // left boundary coefficients
  {0.0, 0.0, 0.0, 0.0, 0.0, -1.5655577299412915, 2.006698028441886, -0.17969802499047988, -0.34551679626695936, 0.031124203444273024, 0.05295031931257198},
  {0.0, 0.0, 0.0, 0.0, -0.4648100847546831, 0.0, 0.31822112431530714, 0.19343747460090951, -0.01748789837432497, -0.029360615787208586, 0.0},
  {0.0, 0.0, 0.0, 0.08895922679561, -0.680116299599615, 0.0, 0.6911055296415899, -0.10654199486276988, 0.006593538025184958, 0.0, 0.0},
  {0.0, 0.0, 0.09023462498760625, -0.21809773186831677, -0.3645866195553873, 0.0, 0.5149776612013978, -0.0225279347652999, 0.0, 0.0, 0.0},
  {0.0, -0.010766469472750029, 0.026116752329614308, 0.07444717053592571, -0.6821178348068156, 0.0, 0.6825798774651727, -0.09025949605114705, 0.0, 0.0, 0.0},
  {-0.016735517732012817, 0.040062884094222144, -0.004209605405399509, 0.02726389466132439, -0.6236615676516715, 0.0, 0.6597484708954711, -0.08246855886193388, 0.0, 0.0, 0.0},
  // central coefficients
  {0.0, 0.0, 0.0, 0.08333333333333333, -0.6666666666666666, 0.0, 0.6666666666666666, -0.08333333333333333, 0.0, 0.0, 0.0},
  // right boundary coefficients
  {0.0, 0.0, 0.0, 0.08246855886193388, -0.6597484708954711, 0.0, 0.6236615676516715, -0.02726389466132439, 0.004209605405399509,-0.040062884094222144, 0.016735517732012817},
  {0.0, 0.0, 0.0, 0.09025949605114705, -0.6825798774651727, 0.0, 0.6821178348068156, -0.07444717053592571, -0.026116752329614308, 0.010766469472750029, 0.0},
  {0.0, 0.0, 0.0, 0.0225279347652999, -0.5149776612013978, 0.0, 0.3645866195553873, 0.21809773186831677, -0.09023462498760625, 0.0, 0.0},
  {0.0, 0.0, -0.006593538025184958, 0.10654199486276988, -0.6911055296415899, 0.0, 0.680116299599615, -0.08895922679561, 0.0, 0.0, 0.0},
  {0.0, 0.029360615787208586, 0.01748789837432497, -0.19343747460090951, -0.31822112431530714, 0.0, 0.4648100847546831, 0.0, 0.0, 0.0, 0.0},
  {-0.05295031931257198, -0.031124203444273024, 0.34551679626695936, 0.17969802499047988, -2.006698028441886, 1.5655577299412915, 0.0, 0.0, 0.0, 0.0, 0.0}
};

/*
Coefficients of the inverse mass/norm matrix.

The central coefficient corresponds to the inner nodes of the grid and is 1.
*/
REAL constant M_INV_B[2*NUM_BOUNDS_B+1] = {
  // left boundary coefficients
  3.131115459882583,
  0.725258121380005,
  1.55005382131324,
  0.817717206132879,
  1.0831139526137645,
  0.9896227063432067,
  // central coefficients
  1.0,
  // right boundary coefficients
  0.9896227063432067,
  1.0831139526137645,
  0.817717206132879,
  1.55005382131324,
  0.725258121380005,
  3.131115459882583
};


