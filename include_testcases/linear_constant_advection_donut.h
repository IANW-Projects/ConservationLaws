// This project is licensed under the terms of the Creative Commons CC BY-NC-ND 4.0 license.

// velocity coefficients a_i
#define a_x 1.0
#define a_y 0.0
#define a_z 0.0

/*
Analytical solution. It is used to calculate the numerical error and related quantities.
*/
inline REAL u_analytical(uint ix, uint iy, uint iz, REAL time) {

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

  return exp(-20*((xx+0.5)*(xx+0.5) + yy*yy + zz*zz));
}

/*
Initial condition.
*/
inline REAL u_init(uint ix, uint iy, uint iz) {

  return u_analytical(ix, iy, iz, (REAL)0);
}

/*
Boundary condition.
*/
inline REAL u_boundary(uint ix, uint iy, uint iz, REAL time) {

  return u_analytical(ix, iy, iz, time);
}
