// This project is licensed under the terms of the Creative Commons CC BY-NC-ND 4.0 license.

//--------------------------------------------------------------------------------------------------
// Functions for field initialisation
//--------------------------------------------------------------------------------------------------

/*
Alfven wave testcase used in [Fluxo](https://github.com/project-fluxo/fluxo), domain [-1,1]^3.
*/

// paramters of the Alfven wave test case
#define CONST_FREQUENCY 1
REAL constant CONST_omega = 6.283185307179586 * CONST_FREQUENCY;
REAL constant CONST_r = 2; // lenght-variable = lenght of computational domain
REAL constant CONST_e = 0.2; // epsilon = 0.2

/*
Analytical solution of the velocity field.
*/
inline REAL4 u_analytical(uint ix, uint iy, uint iz, REAL time) {

  REAL x = (REAL)XMIN + ix*(REAL)DX;
  REAL y = (REAL)YMIN + iy*(REAL)DY;
  REAL z = (REAL)ZMIN + iz*(REAL)DZ;

  REAL nx = rsqrt(CONST_r*CONST_r + 1);
  REAL ny = CONST_r * rsqrt(CONST_r*CONST_r + 1);
  REAL sqr = 1;
  REAL Va = CONST_omega / (ny * sqr);
  REAL phi_alv = CONST_omega / ny * (nx * (x - 0.5*CONST_r) + ny * (y - 0.5*CONST_r)) - Va * time;

  REAL u1 = - CONST_e * ny * cos(phi_alv);
  REAL u2 =   CONST_e * nx * cos(phi_alv);
  REAL u3 =   CONST_e * sin(phi_alv);

  return (REAL4) {u1, u2, u3, (REAL)(0)};
}

/*
Analytical solution of the magnetic field.
*/
inline REAL4 b_analytical(uint ix, uint iy, uint iz, REAL time) {

  REAL x = (REAL)XMIN + ix*(REAL)DX;
  REAL y = (REAL)YMIN + iy*(REAL)DY;
  REAL z = (REAL)ZMIN + iz*(REAL)DZ;

  REAL nx = rsqrt(CONST_r*CONST_r + 1);
  REAL ny = CONST_r * rsqrt(CONST_r*CONST_r + 1);
  REAL sqr = 1;
  REAL Va = CONST_omega / (ny * sqr);
  REAL phi_alv = CONST_omega / ny * (nx * (x - 0.5*CONST_r) + ny * (y - 0.5*CONST_r)) - Va * time;

  REAL B1 = nx + CONST_e * ny * cos(phi_alv) * sqr;
  REAL B2 = ny - CONST_e * nx * cos(phi_alv) * sqr;
  REAL B3 =    - CONST_e * sin(phi_alv) * sqr;

  return (REAL4) {B1, B2, B3, (REAL)(0)};
}

/*
Boundary condition of the magnetic field.
*/
inline REAL4 b_boundary(uint ix, uint iy, uint iz, REAL time) {

  return b_analytical(ix, iy, iz, time);
}

/*
Boundary condition of the velocity field.
*/
inline REAL4 u_boundary(uint ix, uint iy, uint iz, REAL time) {

  return u_analytical(ix, iy, iz, time);
}


/* Boundary condition of the pressure and density */

inline REAL p_boundary(uint ix, uint iy, uint iz, REAL time) {

  return (REAL)(1);
}

inline REAL rho_boundary(uint ix, uint iy, uint iz, REAL time) {

  return (REAL)(1);
}


/* initial conditions for the simulation */

inline REAL4 b_initial(uint ix, uint iy, uint iz, REAL time){

  return b_analytical(ix, iy, iz, time);
}

inline REAL4 u_initial(uint ix, uint iy, uint iz, REAL time){

  return u_analytical(ix, iy, iz, time);
}

inline REAL p_initial(uint idx, uint iy, uint iz, REAL time){
  return (REAL)(1);
}

inline REAL rho_initial(uint idx, uint iy, uint iz, REAL time){

  return (REAL)(1);
}
