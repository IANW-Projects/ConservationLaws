// This project is licensed under the terms of the Creative Commons CC BY-NC-ND 4.0 license.

#ifndef IDEAL_GAS_EULER
#define IDEAL_GAS_EULER


//--------------------------------------------------------------------------------------------------
// Auxiliary functions
//--------------------------------------------------------------------------------------------------

// p = (gamma-1) * (E - 0.5*rho*u^2)
inline REAL compute_energy(REAL p, REAL rho, REAL ux, REAL uy, REAL uz) {

  REAL E = p/(GAMMA-1) + (REAL)(0.5)*rho*(ux*ux+uy*uy+uz*uz);

  return E;
}

// p = (gamma-1) * (E - 0.5*rho*u^2 - 0.5*B^2)
inline REAL compute_pressure(REAL rho, REAL ux, REAL uy, REAL uz, REAL E) {

  REAL p = (GAMMA-1) * (E - (REAL)(0.5)*rho*(ux*ux+uy*uy+uz*uz));

  return p;
}


//--------------------------------------------------------------------------------------------------
// Extended numerical fluxes for the volume terms
//--------------------------------------------------------------------------------------------------


// conservative part
#ifdef USE_FLUX_CENTRAL

  inline void compute_ext_num_flux_x(REAL const* uk, REAL const* um, REAL* ext_num_flux) {

      REAL k_rho = uk[Field_rho];
      REAL k_E   = uk[Field_E];
      REAL k_ux  = uk[Field_ux];
      REAL k_uy  = uk[Field_uy];
      REAL k_uz  = uk[Field_uz];
      REAL k_p   = uk[Field_p];

      ext_num_flux[Field_rho]    = (REAL)(-0.5) * k_rho * k_ux;
      ext_num_flux[Field_rho_ux] = (REAL)(-0.5) * (k_rho*k_ux*k_ux + k_p);
      ext_num_flux[Field_rho_uy] = (REAL)(-0.5) * (k_rho*k_ux*k_uy);
      ext_num_flux[Field_rho_uz] = (REAL)(-0.5) * (k_rho*k_ux*k_uz);
      ext_num_flux[Field_E]      = (REAL)(-0.5) * k_ux * (k_E+k_p);
  }

  inline void compute_ext_num_flux_y(REAL const* uk, REAL const* um, REAL* ext_num_flux) {

      REAL k_rho = uk[Field_rho];
      REAL k_E   = uk[Field_E];
      REAL k_ux  = uk[Field_ux];
      REAL k_uy  = uk[Field_uy];
      REAL k_uz  = uk[Field_uz];
      REAL k_p   = uk[Field_p];

      ext_num_flux[Field_rho]    = (REAL)(-0.5) * k_rho * k_uy;
      ext_num_flux[Field_rho_ux] = (REAL)(-0.5) * (k_rho*k_uy*k_ux);
      ext_num_flux[Field_rho_uy] = (REAL)(-0.5) * (k_rho*k_uy*k_uy + k_p);
      ext_num_flux[Field_rho_uz] = (REAL)(-0.5) * (k_rho*k_uy*k_uz);
      ext_num_flux[Field_E]      = (REAL)(-0.5) * k_uy * (k_E+k_p);
  }

  inline void compute_ext_num_flux_z(REAL const* uk, REAL const* um, REAL* ext_num_flux) {

      REAL k_rho = uk[Field_rho];
      REAL k_E   = uk[Field_E];
      REAL k_ux  = uk[Field_ux];
      REAL k_uy  = uk[Field_uy];
      REAL k_uz  = uk[Field_uz];
      REAL k_p   = uk[Field_p];

      ext_num_flux[Field_rho]    = (REAL)(-0.5) * k_rho * k_uz;
      ext_num_flux[Field_rho_ux] = (REAL)(-0.5) * (k_rho*k_uz*k_ux);
      ext_num_flux[Field_rho_uy] = (REAL)(-0.5) * (k_rho*k_uz*k_uy);
      ext_num_flux[Field_rho_uz] = (REAL)(-0.5) * (k_rho*k_uz*k_uz + k_p);
      ext_num_flux[Field_E]      = (REAL)(-0.5) * k_uz * (k_E+k_p);
  }

#elif defined USE_FLUX_DucrosEtAl

  inline void compute_ext_num_flux_x(REAL const* uk, REAL const* um, REAL* ext_num_flux) {

      REAL m_rho = um[Field_rho];
      REAL m_E   = um[Field_E];
      REAL m_ux  = um[Field_ux];
      REAL m_uy  = um[Field_uy];
      REAL m_uz  = um[Field_uz];
      REAL m_p   = um[Field_p];

      REAL k_rho = uk[Field_rho];
      REAL k_E   = uk[Field_E];
      REAL k_ux  = uk[Field_ux];
      REAL k_uy  = uk[Field_uy];
      REAL k_uz  = uk[Field_uz];
      REAL k_p   = uk[Field_p];

      ext_num_flux[Field_rho]    = (REAL)(-0.25) * (k_rho + m_rho) * (k_ux + m_ux);
      ext_num_flux[Field_rho_ux] = (REAL)(-0.25) * (k_rho*k_ux + m_rho*m_ux) * (k_ux + m_ux) + (REAL)(-0.5) * (k_p + m_p);
      ext_num_flux[Field_rho_uy] = (REAL)(-0.25) * (k_rho*k_uy + m_rho*m_uy) * (k_ux + m_ux);
      ext_num_flux[Field_rho_uz] = (REAL)(-0.25) * (k_rho*k_uz + m_rho*m_uz) * (k_ux + m_ux);
      ext_num_flux[Field_E]      = (REAL)(-0.25) * (k_E + k_p + m_E + m_p) * (k_ux + m_ux);
  }

  inline void compute_ext_num_flux_y(REAL const* uk, REAL const* um, REAL* ext_num_flux) {

      REAL m_rho = um[Field_rho];
      REAL m_E   = um[Field_E];
      REAL m_ux  = um[Field_ux];
      REAL m_uy  = um[Field_uy];
      REAL m_uz  = um[Field_uz];
      REAL m_p   = um[Field_p];

      REAL k_rho = uk[Field_rho];
      REAL k_E   = uk[Field_E];
      REAL k_ux  = uk[Field_ux];
      REAL k_uy  = uk[Field_uy];
      REAL k_uz  = uk[Field_uz];
      REAL k_p   = uk[Field_p];

      ext_num_flux[Field_rho]    = (REAL)(-0.25) * (k_rho + m_rho) * (k_uy + m_uy);
      ext_num_flux[Field_rho_ux] = (REAL)(-0.25) * (k_rho*k_ux + m_rho*m_ux) * (k_uy + m_uy);
      ext_num_flux[Field_rho_uy] = (REAL)(-0.25) * (k_rho*k_uy + m_rho*m_uy) * (k_uy + m_uy) + (REAL)(-0.5) * (k_p + m_p);
      ext_num_flux[Field_rho_uz] = (REAL)(-0.25) * (k_rho*k_uz + m_rho*m_uz) * (k_uy + m_uy);
      ext_num_flux[Field_E]      = (REAL)(-0.25) * (k_E + k_p + m_E + m_p) * (k_uy + m_uy);
  }

  inline void compute_ext_num_flux_z(REAL const* uk, REAL const* um, REAL* ext_num_flux) {

      REAL m_rho = um[Field_rho];
      REAL m_E   = um[Field_E];
      REAL m_ux  = um[Field_ux];
      REAL m_uy  = um[Field_uy];
      REAL m_uz  = um[Field_uz];
      REAL m_p   = um[Field_p];

      REAL k_rho = uk[Field_rho];
      REAL k_E   = uk[Field_E];
      REAL k_ux  = uk[Field_ux];
      REAL k_uy  = uk[Field_uy];
      REAL k_uz  = uk[Field_uz];
      REAL k_p   = uk[Field_p];

      ext_num_flux[Field_rho]    = (REAL)(-0.25) * (k_rho + m_rho) * (k_uz + m_uz);
      ext_num_flux[Field_rho_ux] = (REAL)(-0.25) * (k_rho*k_ux + m_rho*m_ux) * (k_uz + m_uz);
      ext_num_flux[Field_rho_uy] = (REAL)(-0.25) * (k_rho*k_uy + m_rho*m_uy) * (k_uz + m_uz);
      ext_num_flux[Field_rho_uz] = (REAL)(-0.25) * (k_rho*k_uz + m_rho*m_uz) * (k_uz + m_uz) + (REAL)(-0.5) * (k_p + m_p);
      ext_num_flux[Field_E]      = (REAL)(-0.25) * (k_E + k_p + m_E + m_p) * (k_uz + m_uz);
  }

#else

  #error "Error in ideal_MHD.cl: No conservative discretisation specified!"

#endif


//--------------------------------------------------------------------------------------------------
// Field initialisation
//--------------------------------------------------------------------------------------------------

inline void init_fields(uint ix, uint iy, uint iz, global REAL* u) {

  REAL rho = rho_init(ix, iy, iz);
  REAL4 uinit = u_init(ix, iy, iz);
  REAL p = p_init(ix, iy, iz);

  REAL um[NUM_TOTAL_VARS] = {0};

  um[Field_rho] = rho;
  um[Field_rho_ux] = rho*uinit.x;
  um[Field_rho_uy] = rho*uinit.y;
  um[Field_rho_uz] = rho*uinit.z;
  um[Field_E] = compute_energy(p, rho, uinit.x, uinit.y, uinit.z);

  for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
    set_field_component(ix, iy, iz, i, u, um[i]);
  }
}


//--------------------------------------------------------------------------------------------------
// Surface terms
//--------------------------------------------------------------------------------------------------


inline void add_surface_terms(REAL time, uint ix, uint iy, uint iz, global REAL *u, REAL *du_dt) {

  // For periodic boundary conditions and a single block, no surface term has to be used.
#ifndef USE_PERIODIC
  #error "Error in ideal_gas_Euler.cl: Surface terms not implemented for nonperiodic boundaries!"
#endif // USE_PERIODIC

}

//--------------------------------------------------------------------------------------------------
// Computation of auxiliary variables
//---------------------------------------------------------------------------------------------------

inline void compute_auxiliary_variables(REAL time, uint ix, uint iy, uint iz, global REAL *u) {

  REAL um[NUM_TOTAL_VARS] = {0};
  get_field(ix, iy, iz, 0, 0, 0, u, um);

  REAL rho    = um[Field_rho];
  REAL rho_ux = um[Field_rho_ux];
  REAL rho_uy = um[Field_rho_uy];
  REAL rho_uz = um[Field_rho_uz];
  REAL E      = um[Field_E];

  // velocity
  REAL ux = rho_ux / rho;
  REAL uy = rho_uy / rho;
  REAL uz = rho_uz / rho;
  set_field_component(ix, iy, iz, Field_ux, u, ux);
  set_field_component(ix, iy, iz, Field_uy, u, uy);
  set_field_component(ix, iy, iz, Field_uz, u, uz);

  // pressure
  REAL p = compute_pressure(rho, ux, uy, uz, E);
  set_field_component(ix, iy, iz, Field_p, u, p);
}

//--------------------------------------------------------------------------------------------------
// Computation of analytical solution for all conserved variables
//---------------------------------------------------------------------------------------------------

inline void analytical_solution(uint ix, uint iy, uint iz, global REAL *u, REAL time) {

  // No analytical solution
  for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
    set_field_component(ix, iy, iz, i, u, (REAL)(0));
  }
}

#endif // IDEAL_GAS_EULER
