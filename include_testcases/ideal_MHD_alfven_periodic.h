// This project is licensed under the terms of the Creative Commons CC BY-NC-ND 4.0 license.

//--------------------------------------------------------------------------------------------------
// Functions for field initialisation
//--------------------------------------------------------------------------------------------------

/*
Alfven wave testcase used in section 8.6 of
Skew-Symmetric Splitting for Multiscale Gas Dynamics and MHD Turbulence Flows.
Björn Sjögreen, Helen C. Yee, Dmitry Kotov, 2017.
*/

// paramters of the Alfven wave test case
REAL constant CONST_A = 0.1;
REAL constant CONST_alpha = 0.5235987755982988; // pi/6

/*
Analytical solution of the velocity field.
*/
inline REAL4 u_analytical(uint ix, uint iy, uint iz, REAL time) {

  REAL x = (REAL)XMIN + ix*(REAL)DX;
  REAL y = (REAL)YMIN + iy*(REAL)DY;
  REAL z = (REAL)ZMIN + iz*(REAL)DZ;

  REAL cos_alpha = cos(CONST_alpha);
  REAL sin_alpha = sin(CONST_alpha);

  REAL phi = (REAL)(6.283185307179586) * (x*cos_alpha + y*sin_alpha + time);
  REAL cos_phi = cos(phi);
  REAL sin_phi = sin(phi);

  REAL u1 = - CONST_A * sin_alpha * sin_phi;
  REAL u2 = CONST_A * cos_alpha * sin_phi;
  REAL u3 = CONST_A * cos_phi;

  return (REAL4) {u1, u2, u3, (REAL)(0)};
}

/*
Analytical solution of the magnetic field.
*/
inline REAL4 b_analytical(uint ix, uint iy, uint iz, REAL time) {

  REAL x = (REAL)XMIN + ix*(REAL)DX;
  REAL y = (REAL)YMIN + iy*(REAL)DY;
  REAL z = (REAL)ZMIN + iz*(REAL)DZ;

  REAL cos_alpha = cos(CONST_alpha);
  REAL sin_alpha = sin(CONST_alpha);

  REAL phi = (REAL)(6.283185307179586) * (x*cos_alpha + y*sin_alpha + time);
  REAL cos_phi = cos(phi);
  REAL sin_phi = sin(phi);

  REAL B1 = cos_alpha - CONST_A * sin_alpha * sin_phi;
  REAL B2 = sin_alpha + CONST_A * cos_alpha * sin_phi;
  REAL B3 = CONST_A * cos_phi;

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
