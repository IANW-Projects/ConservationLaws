// This project is licensed under the terms of the Creative Commons CC BY-NC-ND 4.0 license.

//--------------------------------------------------------------------------------------------------
// Functions for field initialisation
//--------------------------------------------------------------------------------------------------

/*
Taylor Green vortextestcase used in section 4.2 of
Split form nodal discontinuous Galerkin schemes with summation-by-parts property
for the compressible Euler equations.

Gregor J. Gassner, Andrew R. Winters, David A. Kopriva, 2016.
*/

/*
Initial condition of the density.
*/
inline REAL rho_init(uint ix, uint iy, uint iz) {

	return (REAL)(1);
}

/*
Initial condition of the velocity.
*/
inline REAL4 u_init(uint ix, uint iy, uint iz) {

	REAL x = (REAL)XMIN + ix*(REAL)DX;
	REAL y = (REAL)YMIN + iy*(REAL)DY;
	REAL z = (REAL)ZMIN + iz*(REAL)DZ;

  REAL u1 = sin(x) * cos(y) * cos(z);
  REAL u2 = - cos(x) * sin(y) * cos(z);
  REAL u3 = (REAL)(0);

	return (REAL4) {u1, u2, u3, (REAL)(0)};
}

/*
Initial condition of the pressure.
*/
inline REAL p_init(uint ix, uint iy, uint iz) {

	REAL x = (REAL)XMIN + ix*(REAL)DX;
	REAL y = (REAL)YMIN + iy*(REAL)DY;
	REAL z = (REAL)ZMIN + iz*(REAL)DZ;

	return (REAL)(100.0/GAMMA + 0.0625 * (cos(2*x)*cos(2*z) + 2*cos(2*y) + 2*cos(2*x) + cos(2*y)*cos(2*z)));
}


/*
There is no analytical solution of the density.
*/
inline REAL rho_analytical(uint ix, uint iy, uint iz, REAL time) {

	return rho_init(ix, iy, iz);
}

/*
There is no analytical solution of the velocity.
*/
inline REAL4 u_analytical(uint ix, uint iy, uint iz, REAL time) {

	return u_init(ix, iy, iz);
}

/*
There is no analytical solution of the pressure.
*/
inline REAL p_analytical(uint ix, uint iy, uint iz, REAL time) {

	return p_init(ix, iy, iz);
}
