// This project is licensed under the terms of the Creative Commons CC BY-NC-ND 4.0 license.

#ifndef LINEAR_VARIABLE_ADVECTION_H
#define LINEAR_VARIABLE_ADVECTION_H

/*
The linear advection equation
  \partial_t u + div(a u) = 0
with variable velocity $a = (a_x, a_y, a_z)$.
*/


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

inline void add_surface_terms(REAL time, uint ix, uint iy, uint iz, global REAL const *u, REAL *du_dt) {

  // For periodic boundary conditions and a single block, no surface term has to be used.
  #ifdef NON_PERIODIC_BOUNDARY_EXISTS

    REAL um[NUM_TOTAL_VARS] = {0};
    get_field(ix, iy, iz, 0, 0, 0, u, um);

    REAL u_bound = u_boundary(ix, iy, iz, time);

    // Apply the usual upwind numerical flux at the boundaries.
  #ifndef USE_PERIODIC_X
    du_dt[Field_u] += + (REAL)(M_INV_X[0]/DX) *
                        ((check_bound_l(ix,1) * (um[Field_ax] > 0)) * (REAL)(-1) + (check_bound_xr(ix,1) * (um[Field_ax] < 0)) * (REAL)(1)) * um[Field_ax] * (um[Field_u] - u_bound)
  #endif
  #ifdef USE_PERIODIC_Y
    du_dt[Field_u] += + (REAL)(M_INV_Y[0]/DY) *
                        ((check_bound_l(iy,1) * (um[Field_ay] > 0)) * (REAL)(-1) + (check_bound_yr(iy,1) * (um[Field_ay] < 0)) * (REAL)(1)) * um[Field_ay] * (um[Field_u] - u_bound)
  #endif
  #ifdef USE_PERIODIC_Z
    du_dt[Field_u] += + (REAL)(M_INV_Z[0]/DZ) *
                        ((check_bound_l(iz,1) * (um[Field_az] > 0)) * (REAL)(-1) + (check_bound_zr(iz,1) * (um[Field_az] < 0)) * (REAL)(1)) * um[Field_az] * (um[Field_u] - u_bound);
  #endif
  #endif // NON_PERIODIC_BOUNDARY_EXISTS

  return;
}

//--------------------------------------------------------------------------------------------------
// Computation of auxiliary variables
//---------------------------------------------------------------------------------------------------

inline void compute_auxiliary_variables(REAL time, uint ix, uint iy, uint iz, global REAL *u) {

  REAL x = (REAL)XMIN + ix*(REAL)DX;
  REAL y = (REAL)YMIN + iy*(REAL)DY;

  set_field_component(ix, iy, iz, Field_ax, u, -y);
  set_field_component(ix, iy, iz, Field_ay, u, x);
  set_field_component(ix, iy, iz, Field_az, u, (REAL)(0));
}

//--------------------------------------------------------------------------------------------------
// Computation of analytical solution for all conserved variables
//---------------------------------------------------------------------------------------------------

inline void analytical_solution(uint ix, uint iy, uint iz, global REAL *u, REAL time) {

  REAL u_ana =  u_analytical(ix, iy, iz, time);

  set_field_component(ix, iy, iz, Field_u, u, u_ana);
}

#endif // LINEAR_VARIABLE_ADVECTION_H
