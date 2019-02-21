// This project is licensed under the terms of the Creative Commons CC BY-NC-ND 4.0 license.

//--------------------------------------------------------------------------------------------------
// Functions for field initialisation
//--------------------------------------------------------------------------------------------------

/*
Constant velocity and pressure result in the linear advection equation for the density.
*/

// velocity coefficients a_i
#define a_x M_PI
#define a_y M_PI
#define a_z M_PI


/*
Analytical solution of the density.
*/
inline REAL rho_analytical(uint ix, uint iy, uint iz, REAL time) {

 // Coordinates of (ix,iy,iz)
  REAL x = (REAL)XMIN + ix*(REAL)DX;
  REAL y = (REAL)YMIN + iy*(REAL)DY;
  REAL z = (REAL)ZMIN + iz*(REAL)DZ;

  // The solution is u(t,x) = u_0(x - a*t).
  REAL dx = -time*a_x;
  REAL dy = -time*a_y;
  REAL dz = -time*a_z;

  // A periodic solution is desired.
  // Note: fmod(x, y) seems to return a value in [-y,y] (depending on the sign of x) but we want to have a value in [0,y].
  REAL xx = fmod(x + dx - (REAL)XMIN, (REAL)(XMAX - XMIN)); xx = xx + (xx < 0) * (REAL)(XMAX - XMIN) + (REAL)XMIN;
  REAL yy = fmod(y + dy - (REAL)YMIN, (REAL)(YMAX - YMIN)); yy = yy + (yy < 0) * (REAL)(YMAX - YMIN) + (REAL)YMIN;
  REAL zz = fmod(z + dz - (REAL)ZMIN, (REAL)(ZMAX - ZMIN)); zz = zz + (zz < 0) * (REAL)(ZMAX - ZMIN) + (REAL)ZMIN;

	return 1 + sin(0.5*xx) * sin(0.5*yy) * sin(0.5*zz);
}

/*
Analytical solution of the velocity.
*/
inline REAL4 u_analytical(uint ix, uint iy, uint iz, REAL time) {

	return (REAL4) {a_x, a_y, a_z, 0};
}

/*
Analytical solution of the pressure.
*/
inline REAL p_analytical(uint ix, uint iy, uint iz, REAL time) {

	return (REAL)(10);
}


/*
Boundary condition of the density.
*/
inline REAL rho_boundary(uint ix, uint iy, uint iz, REAL time) {

	return rho_analytical(ix, iy, iz, time);
}

/*
Boundary condition of the velocity.
*/
inline REAL4 u_boundary(uint ix, uint iy, uint iz, REAL time) {

	return u_analytical(ix, iy, iz, time);
}

/*
Boundary condition of the pressure.
*/
inline REAL p_boundary(uint ix, uint iy, uint iz, REAL time) {

	return p_analytical(ix, iy, iz, time);
}


/*
Initial condition of the density.
*/
inline REAL rho_init(uint ix, uint iy, uint iz) {

	return rho_analytical(ix, iy, iz, (REAL)(0));
}

/*
Initial condition of the velocity.
*/
inline REAL4 u_init(uint ix, uint iy, uint iz) {

	return u_analytical(ix, iy, iz, (REAL)(0));
}

/*
Initial condition of the pressure.
*/
inline REAL p_init(uint ix, uint iy, uint iz) {

	return p_analytical(ix, iy, iz, (REAL)(0));
}
