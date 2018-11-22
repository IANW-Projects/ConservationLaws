// This project is licensed under the terms of the Creative Commons CC BY-NC-ND 4.0 license.

#ifndef LINEAR_CONSTANT_ADVECTION_H
#define LINEAR_CONSTANT_ADVECTION_H

/*
The linear advection equation
  \partial_t u + div(a u) = 0
with constant velocity $a = (a_x, a_y, a_z)$.
*/

// velocity coefficients a_i
#define a_x 1.0
#define a_y 0.0
#define a_z 0.0


inline void compute_auxiliary_variables(REAL time, uint ix, uint iy, uint iz, global REAL *u) {

  // nothing to do
  return;
}


inline void add_surface_terms(REAL time, uint ix, uint iy, uint iz, global REAL *u, REAL *du_dt) {

  // nothing to do in a periodic domain
  return;
}


// (extended) numerical fluxes for the volume terms
inline void compute_ext_num_flux_x(REAL const* uk, REAL const* um, REAL* ext_num_flux) {
  ext_num_flux[Field_u] = (REAL)(-0.5) * a_x * (uk[Field_u] + um[Field_u]);
}


inline void compute_ext_num_flux_y(REAL const* uk, REAL const* um, REAL* ext_num_flux) {
  ext_num_flux[Field_u] = (REAL)(-0.5) * a_y * (uk[Field_u] + um[Field_u]);
}


inline void compute_ext_num_flux_z(REAL const* uk, REAL const* um, REAL* ext_num_flux) {
  ext_num_flux[Field_u] = (REAL)(-0.5) * a_z * (uk[Field_u] + um[Field_u]);
}


//TODO: (extended) numerical fluxes for the surface terms: upwind

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
Initial condition.
*/
inline REAL u_init(uint ix, uint iy, uint iz) {

	return u_analytical(ix, iy, iz, (REAL)0);
}

// Initialise fields.
inline void init_fields(uint ix, uint iy, uint iz, global REAL* u) {

  REAL uinit = u_init(ix, iy, iz);

  set_field_component(ix, iy, iz, Field_u, u, uinit);
}


#endif // LINEAR_CONSTANT_ADVECTION_H
