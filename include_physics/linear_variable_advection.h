// This project is licensed under the terms of the Creative Commons CC BY-NC-ND 4.0 license.

#ifndef LINEAR_VARIABLE_ADVECTION_H
#define LINEAR_VARIABLE_ADVECTION_H

/*
The linear advection equation
  \partial_t u + div(a u) = 0
with variable velocity $a = (a_x, a_y, a_z)$.
*/


inline void compute_auxiliary_variables(REAL time, uint ix, uint iy, uint iz, global REAL *u) {

  REAL x = (REAL)XMIN + ix*(REAL)DX;
	REAL y = (REAL)YMIN + iy*(REAL)DY;

  set_field_component(ix, iy, iz, Field_ax, u, -y);
  set_field_component(ix, iy, iz, Field_ay, u, x);
  set_field_component(ix, iy, iz, Field_az, u, (REAL)(0));
}


inline void add_surface_terms(REAL time, uint ix, uint iy, uint iz, global REAL const *u, REAL *du_dt) {

  // nothing to do in a periodic domain
  return;
}


// extended numerical fluxes for the volume terms
inline void compute_ext_num_flux_x(REAL const* uk, REAL const* um, REAL* ext_num_flux) {

  #ifdef USE_CENTRAL_FLUX
    ext_num_flux[Field_u] = (REAL)(-0.5) * (uk[Field_ax]*uk[Field_u] + um[Field_ax]*um[Field_u]);
	#elif defined USE_SPLIT_FLUX
    ext_num_flux[Field_u] = (REAL)(-0.25) * (uk[Field_ax] + um[Field_ax]) * (uk[Field_u] + um[Field_u]);
	#elif defined USE_PRODUCT_FLUX
    ext_num_flux[Field_u] = (REAL)(-0.5) * (uk[Field_ax]*um[Field_u] + um[Field_ax]*uk[Field_u]);
	#else
    #error "Error in linear_variable_advection.cl: No discretisation specified!"
  #endif
}

inline void compute_ext_num_flux_y(REAL const* uk, REAL const* um, REAL* ext_num_flux) {

  #ifdef USE_CENTRAL_FLUX
    ext_num_flux[Field_u] = (REAL)(-0.5) * (uk[Field_ay]*uk[Field_u] + um[Field_ay]*um[Field_u]);
	#elif defined USE_SPLIT_FLUX
    ext_num_flux[Field_u] = (REAL)(-0.25) * (uk[Field_ay] + um[Field_ay]) * (uk[Field_u] + um[Field_u]);
	#elif defined USE_PRODUCT_FLUX
    ext_num_flux[Field_u] = (REAL)(-0.5) * (uk[Field_ay]*um[Field_u] + um[Field_ay]*uk[Field_u]);
	#else
    #error "Error in linear_variable_advection.cl: No discretisation specified!"
  #endif
}

inline void compute_ext_num_flux_z(REAL const* uk, REAL const* um, REAL* ext_num_flux) {

  #ifdef USE_CENTRAL_FLUX
    ext_num_flux[Field_u] = (REAL)(-0.5) * (uk[Field_az]*uk[Field_u] + um[Field_az]*um[Field_u]);
	#elif defined USE_SPLIT_FLUX
    ext_num_flux[Field_u] = (REAL)(-0.25) * (uk[Field_az] + um[Field_az]) * (uk[Field_u] + um[Field_u]);
	#elif defined USE_PRODUCT_FLUX
    ext_num_flux[Field_u] = (REAL)(-0.5) * (uk[Field_az]*um[Field_u] + um[Field_az]*uk[Field_u]);
	#else
    #error "Error in linear_variable_advection.cl: No discretisation specified!"
  #endif
}



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

// Initialise fields.
inline void init_fields(uint ix, uint iy, uint iz, global REAL* u) {

  REAL uinit = u_init(ix, iy, iz);

  set_field_component(ix, iy, iz, Field_u, u, uinit);
}

#endif // LINEAR_VARIABLE_ADVECTION_H
