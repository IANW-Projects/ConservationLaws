// This project is licensed under the terms of the Creative Commons CC BY-NC-ND 4.0 license.

/*
Analytical solution. It is used to calculate the numerical error and related quantities.
*/
inline REAL u_analytical(uint ix, uint iy, uint iz, REAL time) {

	REAL x = (REAL)XMIN + ix*(REAL)DX;
	REAL y = (REAL)YMIN + iy*(REAL)DY;
	REAL z = (REAL)ZMIN + iz*(REAL)DZ;

	REAL s = sin(time);
  REAL c = cos(time);

	REAL xx = c*x + s*y;
	REAL yy = -s*x + c*y;
	REAL zz = z;

	REAL fac = exp(-20*((xx-(REAL)(0.5))*(xx-(REAL)(0.5)) + yy*yy));
	REAL Bx0 = -4*yy*fac;
	REAL By0 = (4*xx-2)*fac;

	return (REAL) sqrt((c*Bx0-s*By0)*(c*Bx0-s*By0) + (s*Bx0+c*By0)*(s*Bx0+c*By0));
}

/*
Initial condition of the magnetic field.
*/
inline REAL u_init(uint ix, uint iy, uint iz) {

	return u_analytical(ix, iy, iz, (REAL)0);
}
