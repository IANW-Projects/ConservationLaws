// This project is licensed under the terms of the Creative Commons CC BY-NC-ND 4.0 license.

/* Contains the coefficients of the sixth order diagonal norm derivative operator of
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

#ifdef USE_PERIODIC
//------------------------------------------------------------------------------
// Periodic Derivative Operator
//------------------------------------------------------------------------------
#define ORDER 6
#define NUM_BOUNDS 0
#define STENCIL_WIDTH 7

/*
Coefficients of the first derivative operator.

For each row, the central coefficient corresponds to the nodes at which the
derivative is computed.
*/
REAL constant SBP_diff[2*NUM_BOUNDS+1][STENCIL_WIDTH] = {
  // central coefficients
  {-0.016666666666666666, 0.15, -0.75, 0.0, 0.75, -0.15, 0.016666666666666666},
};


/*
Coefficients of the inverse mass/norm matrix.

The central coefficient corresponds to the inner nodes of the grid and is 1.
*/
REAL constant M_INV[2*NUM_BOUNDS+1] = {
  // central coefficients
  1.0,
};

#else
//------------------------------------------------------------------------------
// Non-Periodic Derivative Operator
//------------------------------------------------------------------------------

#define ORDER 6
#define NUM_BOUNDS 6
#define STENCIL_WIDTH 11

/*
Coefficients of the first derivative operator.

For each row, the central coefficient corresponds to the nodes at which the
derivative is computed.
*/
REAL constant SBP_diff[2*NUM_BOUNDS+1][STENCIL_WIDTH] = {
  // left boundary coefficients
  {0.0, 0.0, 0.0, 0.0, 0.0, -1.5825335189391163, 1.905066305223826, 0.3717366351625272, -1.2202725474393727, 0.6177375631914426, -0.09173443719930642},
  {0.0, 0.0, 0.0, 0.0, -0.4329018563223175, 0.0, -0.0043147701101584396, 0.8419628735536502, -0.5064721551652377, 0.10172590804406338, 0.0},
  {0.0, 0.0, 0.0, -0.1871572605434649, 0.009559818025328907, 0.0, -0.6857863027173244, 1.2691196360506578, -0.40573589081519734, 0.0, 0.0},
  {0.0, 0.0, 0.31079492442619894, -0.9436928531442433, 0.3469241773962804, 0.0, 0.19345960067176712, 0.0790788082353673, 0.013435342414629596, 0.0, 0.0},
  {0.0, -0.21407896407261648, 0.772407007744065, -0.8735770809529855, -0.26323473403580044, 0.0, 0.7247323431086286, -0.1645296432652025, 0.01828107147391139, 0.0, 0.0},
  {0.02858572483124434, -0.1394983371764724, 0.2511244035524303, -0.09675197674330115, -0.6516651065805195, 0.0, 0.7397091390607521, -0.1479418278121504, 0.016437980868016712, 0.0, 0.0},
  // central coefficients
  {0.0, 0.0, -0.016666666666666666, 0.15, -0.75, 0.0, 0.75, -0.15, 0.016666666666666666, 0.0, 0.0},
  // right boundary coefficients
  {0.0, 0.0, -0.016437980868016712, 0.1479418278121504, -0.7397091390607521, 0.0, 0.6516651065805195, 0.09675197674330115,-0.2511244035524303, 0.1394983371764724, -0.02858572483124434},
  {0.0, 0.0, -0.01828107147391139, 0.1645296432652025, -0.7247323431086286, 0.0, 0.26323473403580044, 0.8735770809529855, -0.772407007744065, 0.21407896407261648, 0.0},
  {0.0, 0.0, -0.013435342414629596, -0.0790788082353673, -0.19345960067176712, 0.0, -0.3469241773962804, 0.9436928531442433, -0.31079492442619894, 0.0, 0.0},
  {0.0, 0.0, 0.40573589081519734, -1.2691196360506578, 0.6857863027173244, 0.0, -0.009559818025328907, 0.1871572605434649,0.0, 0.0, 0.0},
  {0.0, -0.10172590804406338, 0.5064721551652377, -0.8419628735536502, 0.0043147701101584396, 0.0, 0.4329018563223175, 0.0, 0.0, 0.0, 0.0},
  {0.09173443719930642, -0.6177375631914426, 1.2202725474393727, -0.3717366351625272, -1.905066305223826, 1.5825335189391163, 0.0, 0.0, 0.0, 0.0, 0.0}
};


/*
Coefficients of the inverse mass/norm matrix.

The central coefficient corresponds to the inner nodes of the grid and is 1.
*/
REAL constant M_INV[2*NUM_BOUNDS+1] = {
  // left boundary coefficients
  3.1650670378782326,
  0.719220844085574,
  1.5935079306528956,
  0.8061205448777757,
  1.0968642884346833,
  0.9862788520810027,
  // central coefficients
  1.0,
  // right boundary coefficients
  0.9862788520810027,
  1.0968642884346833,
  0.8061205448777757,
  1.5935079306528956,
  0.719220844085574,
  3.1650670378782326
};

#endif // USE_PERIODIC
