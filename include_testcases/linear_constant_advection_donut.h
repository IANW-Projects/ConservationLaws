// This project is licensed under the terms of the Creative Commons CC BY-NC-ND 4.0 license.

// velocity coefficients a_i
#define a_x 0.0
#define a_y 0.0
#define a_z 2.0

/*
Analytical solution. It is used to calculate the numerical error and related quantities.
*/
inline REAL u_analytical(uint ix, uint iy, uint iz, REAL time) {

	REAL x = (REAL)XMIN + ix*(REAL)DX;
	REAL y = (REAL)YMIN + iy*(REAL)DY;
	REAL z = (REAL)ZMIN + iz*(REAL)DZ;

	REAL s = 0;
  	REAL c = 1;

	REAL xx = c*x + s*y;
	REAL yy = -s*x + c*y;
	REAL zz = z;

	REAL dx = -(REAL)0.5 + time*a_x;
	REAL dy =  time*a_y;

	REAL fac = exp(-20*((xx+dx)*(xx+dx) + (yy+dy)*(yy+dy)));
	REAL Bx0 = -4*(yy+dy)*fac;
	REAL By0 = (4*(xx+dx))*fac;

	return (REAL) sqrt((c*Bx0-s*By0)*(c*Bx0-s*By0) + (s*Bx0+c*By0)*(s*Bx0+c*By0));
}

/*
Initial condition.
*/
inline REAL u_init(uint ix, uint iy, uint iz) {

	return u_analytical(ix, iy, iz, (REAL)0);
}
