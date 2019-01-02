// This project is licensed under the terms of the Creative Commons CC BY-NC-ND 4.0 license.

#ifndef LINEAR_CONSTANT_ADVECTION_H
#define LINEAR_CONSTANT_ADVECTION_H

/*
The linear advection equation
  \partial_t u + div(a u) = 0
with constant velocity $a = (a_x, a_y, a_z)$.
*/



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


//--------------------------------------------------------------------------------------------------
// Field initialisation
//--------------------------------------------------------------------------------------------------

// Initialise fields.
inline void init_fields(uint ix, uint iy, uint iz, global REAL* u) {

  REAL uinit = u_init(ix, iy, iz);

  set_field_component(ix, iy, iz, Field_u, u, uinit);
}

//--------------------------------------------------------------------------------------------------
// Surface terms
//--------------------------------------------------------------------------------------------------

inline void add_surface_terms(REAL time, uint ix, uint iy, uint iz, global REAL *u, REAL *du_dt) {

  // For periodic boundary conditions and a single block, no surface term has to be used.
  #ifndef USE_PERIODIC

    REAL um[NUM_TOTAL_VARS] = {0};
    get_field(ix, iy, iz, 0, 0, 0, u, um);

    REAL u_bound = u_boundary(ix, iy, iz, time);

    // Apply the usual upwind numerical flux at the boundaries.
    du_dt[Field_u] += + (REAL)(M_INV[0]/DX) *
                        ((check_bound_l(ix,1) * (a_x > 0)) * (REAL)(-1) + (check_bound_xr(ix,1) * (a_x < 0)) * (REAL)(1)) * a_x * (um[Field_u] - u_bound)
                      + (REAL)(M_INV[0]/DY) *
                        ((check_bound_l(iy,1) * (a_y > 0)) * (REAL)(-1) + (check_bound_yr(iy,1) * (a_y < 0)) * (REAL)(1)) * a_y * (um[Field_u] - u_bound)
                      + (REAL)(M_INV[0]/DZ) *
                        ((check_bound_l(iz,1) * (a_z > 0)) * (REAL)(-1) + (check_bound_zr(iz,1) * (a_z < 0)) * (REAL)(1)) * a_z * (um[Field_u] - u_bound);

  #endif // USE_PERIODIC

  return;
}

//--------------------------------------------------------------------------------------------------
// Computation of auxiliary variables
//---------------------------------------------------------------------------------------------------

inline void compute_auxiliary_variables(REAL time, uint ix, uint iy, uint iz, global REAL *u) {

  // nothing to do
  return;
}

//--------------------------------------------------------------------------------------------------
// Computation of analytical solution for all conserved variables
//---------------------------------------------------------------------------------------------------

inline void analytical_solution(uint ix, uint iy, uint iz, global REAL *u, REAL time) {

  REAL u_ana =  u_analytical(ix, iy, iz, time);

  set_field_component(ix, iy, iz, Field_u, u, u_ana);
}


#endif // LINEAR_CONSTANT_ADVECTION_H
