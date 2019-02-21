// This project is licensed under the terms of the Creative Commons CC BY-NC-ND 4.0 license.

#ifndef IDEAL_MHD
#define IDEAL_MHD

// ideal gas constant
REAL constant GAMMA = (REAL)(1.4);


//--------------------------------------------------------------------------------------------------
// Auxiliary functions
//--------------------------------------------------------------------------------------------------

// p = (gamma-1) * (E - 0.5*rho*u^2 - 0.5*B^2)
inline REAL compute_energy(REAL p, REAL rho, REAL ux, REAL uy, REAL uz, REAL Bx, REAL By, REAL Bz) {

  REAL E = p/(GAMMA-1) + (REAL)(0.5)*rho*(ux*ux+uy*uy+uz*uz) + (REAL)(0.5)*(Bx*Bx+By*By+Bz*Bz);

  return E;
}

// p = (gamma-1) * (E - 0.5*rho*u^2 - 0.5*B^2)
inline REAL compute_pressure(REAL rho, REAL ux, REAL uy, REAL uz, REAL E, REAL Bx, REAL By, REAL Bz) {

  REAL p = (GAMMA-1) * (E - (REAL)(0.5)*rho*(ux*ux+uy*uy+uz*uz) - (REAL)(0.5)*(Bx*Bx+By*By+Bz*Bz));

  return p;
}


//--------------------------------------------------------------------------------------------------
// Extended numerical fluxes for the volume terms
//--------------------------------------------------------------------------------------------------


// conservative part
#ifdef USE_FLUX_CENTRAL

  inline void add_ext_num_flux_x_conservative(REAL const* uk, REAL const* um, REAL* ext_num_flux) {

      REAL k_rho = uk[Field_rho];
      REAL k_E   = uk[Field_E];
      REAL k_Bx  = uk[Field_Bx];
      REAL k_By  = uk[Field_By];
      REAL k_Bz  = uk[Field_Bz];
      REAL k_ux  = uk[Field_ux];
      REAL k_uy  = uk[Field_uy];
      REAL k_uz  = uk[Field_uz];
      REAL k_p   = uk[Field_p];

      REAL k_B2 = k_Bx*k_Bx + k_By*k_By + k_Bz*k_Bz;

      ext_num_flux[Field_rho]    += (REAL)(-0.5) * k_rho * k_ux;
      ext_num_flux[Field_rho_ux] += (REAL)(-0.5) * (k_rho*k_ux*k_ux + k_p + (REAL)(0.5)*k_B2 - k_Bx*k_Bx);
      ext_num_flux[Field_rho_uy] += (REAL)(-0.5) * (k_rho*k_ux*k_uy - k_Bx*k_By);
      ext_num_flux[Field_rho_uz] += (REAL)(-0.5) * (k_rho*k_ux*k_uz - k_Bx*k_Bz);
      ext_num_flux[Field_E]      += (REAL)(-0.5) * (k_ux*(k_E+k_p+(REAL)(0.5)*k_B2) - k_Bx*(k_ux*k_Bx + k_uy*k_By + k_uz*k_Bz));
      ext_num_flux[Field_Bx]     += (REAL)(-0.5) * (0);
      ext_num_flux[Field_By]     += (REAL)(-0.5) * (k_ux*k_By - k_uy*k_Bx);
      ext_num_flux[Field_Bz]     += (REAL)(-0.5) * (k_ux*k_Bz - k_uz*k_Bx);
  }

  inline void add_ext_num_flux_y_conservative(REAL const* uk, REAL const* um, REAL* ext_num_flux) {

      REAL k_rho = uk[Field_rho];
      REAL k_E   = uk[Field_E];
      REAL k_Bx  = uk[Field_Bx];
      REAL k_By  = uk[Field_By];
      REAL k_Bz  = uk[Field_Bz];
      REAL k_ux  = uk[Field_ux];
      REAL k_uy  = uk[Field_uy];
      REAL k_uz  = uk[Field_uz];
      REAL k_p   = uk[Field_p];

      REAL k_B2 = k_Bx*k_Bx + k_By*k_By + k_Bz*k_Bz;

      ext_num_flux[Field_rho]    += (REAL)(-0.5) * k_rho * k_uy;
      ext_num_flux[Field_rho_ux] += (REAL)(-0.5) * (k_rho*k_uy*k_ux - k_By*k_Bx);
      ext_num_flux[Field_rho_uy] += (REAL)(-0.5) * (k_rho*k_uy*k_uy + k_p + (REAL)(0.5)*k_B2 - k_By*k_By);
      ext_num_flux[Field_rho_uz] += (REAL)(-0.5) * (k_rho*k_uy*k_uz - k_By*k_Bz);
      ext_num_flux[Field_E]      += (REAL)(-0.5) * (k_uy*(k_E+k_p+(REAL)(0.5)*k_B2) - k_By*(k_ux*k_Bx + k_uy*k_By + k_uz*k_Bz));
      ext_num_flux[Field_Bx]     += (REAL)(-0.5) * (k_uy*k_Bx - k_ux*k_By);
      ext_num_flux[Field_By]     += (REAL)(-0.5) * (0);
      ext_num_flux[Field_Bz]     += (REAL)(-0.5) * (k_uy*k_Bz - k_uz*k_By);
  }

  inline void add_ext_num_flux_z_conservative(REAL const* uk, REAL const* um, REAL* ext_num_flux) {

      REAL k_rho = uk[Field_rho];
      REAL k_E   = uk[Field_E];
      REAL k_Bx  = uk[Field_Bx];
      REAL k_By  = uk[Field_By];
      REAL k_Bz  = uk[Field_Bz];
      REAL k_ux  = uk[Field_ux];
      REAL k_uy  = uk[Field_uy];
      REAL k_uz  = uk[Field_uz];
      REAL k_p   = uk[Field_p];

      REAL k_B2 = k_Bx*k_Bx + k_By*k_By + k_Bz*k_Bz;

      ext_num_flux[Field_rho]    += (REAL)(-0.5) * k_rho * k_uz;
      ext_num_flux[Field_rho_ux] += (REAL)(-0.5) * (k_rho*k_uz*k_ux - k_Bz*k_Bx);
      ext_num_flux[Field_rho_uy] += (REAL)(-0.5) * (k_rho*k_uz*k_uy - k_Bz*k_By);
      ext_num_flux[Field_rho_uz] += (REAL)(-0.5) * (k_rho*k_uz*k_uz + k_p + (REAL)(0.5)*k_B2 - k_Bz*k_Bz);
      ext_num_flux[Field_E]      += (REAL)(-0.5) * (k_uz*(k_E+k_p+(REAL)(0.5)*k_B2) - k_Bz*(k_ux*k_Bx + k_uy*k_By + k_uz*k_Bz));
      ext_num_flux[Field_Bx]     += (REAL)(-0.5) * (k_uz*k_Bx - k_ux*k_Bz);
      ext_num_flux[Field_By]     += (REAL)(-0.5) * (k_uz*k_By - k_uy*k_Bz);
      ext_num_flux[Field_Bz]     += (REAL)(-0.5) * (0);
  }

#elif defined USE_FLUX_SJOGREEN_YEE_DS2

  inline void add_ext_num_flux_x_conservative(REAL const* uk, REAL const* um, REAL* ext_num_flux) {

      REAL m_rho = um[Field_rho];
      REAL m_E   = um[Field_E];
      REAL m_Bx  = um[Field_Bx];
      REAL m_By  = um[Field_By];
      REAL m_Bz  = um[Field_Bz];
      REAL m_ux  = um[Field_ux];
      REAL m_uy  = um[Field_uy];
      REAL m_uz  = um[Field_uz];
      REAL m_p   = um[Field_p];

      REAL m_B2 = m_Bx*m_Bx + m_By*m_By + m_Bz*m_Bz;

      REAL k_rho = uk[Field_rho];
      REAL k_E   = uk[Field_E];
      REAL k_Bx  = uk[Field_Bx];
      REAL k_By  = uk[Field_By];
      REAL k_Bz  = uk[Field_Bz];
      REAL k_ux  = uk[Field_ux];
      REAL k_uy  = uk[Field_uy];
      REAL k_uz  = uk[Field_uz];
      REAL k_p   = uk[Field_p];

      REAL k_B2 = k_Bx*k_Bx + k_By*k_By + k_Bz*k_Bz;

      ext_num_flux[Field_rho]    += (REAL)(-0.25) * (k_rho + m_rho) * (k_ux + m_ux);
      ext_num_flux[Field_rho_ux] += (REAL)(-0.25) * (k_rho*k_ux + m_rho*m_ux) * (k_ux + m_ux)
                                  + (REAL)(-0.5) * (k_p + (REAL)(0.5) * k_B2)
                                  - (REAL)(-0.25) * (k_Bx + m_Bx) * (k_Bx + m_Bx);
      ext_num_flux[Field_rho_uy] += (REAL)(-0.25) * (k_rho*k_uy + m_rho*m_uy) * (k_ux + m_ux)
                                  - (REAL)(-0.25) * (k_Bx + m_Bx) * (k_By + m_By);
      ext_num_flux[Field_rho_uz] += (REAL)(-0.25) * (k_rho*k_uz + m_rho*m_uz) * (k_ux + m_ux)
                                  - (REAL)(-0.25) * (k_Bx + m_Bx) * (k_Bz + m_Bz);
      ext_num_flux[Field_E]      += (REAL)(-0.25) * (k_E + k_p + (REAL)(0.5)*k_B2 + m_E + m_p + (REAL)(0.5)*m_B2) * (k_ux + m_ux)
                                  - (REAL)(-0.5) * (k_Bx*(k_ux*k_Bx + k_uy*k_By + k_uz*k_Bz));
      ext_num_flux[Field_Bx]     += (REAL)(-0.5) * (0);
      ext_num_flux[Field_By]     += (REAL)(-0.25) * (k_ux + m_ux) * (k_By + m_By)
                                  - (REAL)(-0.25) * (k_uy + m_uy) * (k_Bx + m_Bx);
      ext_num_flux[Field_Bz]     += (REAL)(-0.25) * (k_ux + m_ux) * (k_Bz + m_Bz)
                                  - (REAL)(-0.25) * (k_uz + m_uz) * (k_Bx + m_Bx);
  }

  inline void add_ext_num_flux_y_conservative(REAL const* uk, REAL const* um, REAL* ext_num_flux) {

      REAL m_rho = um[Field_rho];
      REAL m_E   = um[Field_E];
      REAL m_Bx  = um[Field_Bx];
      REAL m_By  = um[Field_By];
      REAL m_Bz  = um[Field_Bz];
      REAL m_ux  = um[Field_ux];
      REAL m_uy  = um[Field_uy];
      REAL m_uz  = um[Field_uz];
      REAL m_p   = um[Field_p];

      REAL m_B2 = m_Bx*m_Bx + m_By*m_By + m_Bz*m_Bz;

      REAL k_rho = uk[Field_rho];
      REAL k_E   = uk[Field_E];
      REAL k_Bx  = uk[Field_Bx];
      REAL k_By  = uk[Field_By];
      REAL k_Bz  = uk[Field_Bz];
      REAL k_ux  = uk[Field_ux];
      REAL k_uy  = uk[Field_uy];
      REAL k_uz  = uk[Field_uz];
      REAL k_p   = uk[Field_p];

      REAL k_B2 = k_Bx*k_Bx + k_By*k_By + k_Bz*k_Bz;

      ext_num_flux[Field_rho]    += (REAL)(-0.25) * (k_rho + m_rho) * (k_uy + m_uy);
      ext_num_flux[Field_rho_ux] += (REAL)(-0.25) * (k_rho*k_ux + m_rho*m_ux) * (k_uy + m_uy)
                                  - (REAL)(-0.25) * (k_By + m_By) * (k_Bx + m_Bx);
      ext_num_flux[Field_rho_uy] += (REAL)(-0.25) * (k_rho*k_uy + m_rho*m_uy) * (k_uy + m_uy)
                                  + (REAL)(-0.5) * (k_p + (REAL)(0.5) * k_B2)
                                  - (REAL)(-0.25) * (k_By + m_By) * (k_By + m_By);
      ext_num_flux[Field_rho_uz] += (REAL)(-0.25) * (k_rho*k_uz + m_rho*m_uz) * (k_uy + m_uy)
                                  - (REAL)(-0.25) * (k_By + m_By) * (k_Bz + m_Bz);
      ext_num_flux[Field_E]      += (REAL)(-0.25) * (k_E + k_p + (REAL)(0.5)*k_B2 + m_E + m_p + (REAL)(0.5)*m_B2) * (k_uy + m_uy)
                                  - (REAL)(-0.5) * (k_By*(k_ux*k_Bx + k_uy*k_By + k_uz*k_Bz));
      ext_num_flux[Field_Bx]     += (REAL)(-0.25) * (k_uy + m_uy) * (k_Bx + m_Bx)
                                  - (REAL)(-0.25) * (k_ux + m_ux) * (k_By + m_By);
      ext_num_flux[Field_By]     += (REAL)(-0.5) * (0);
      ext_num_flux[Field_Bz]     += (REAL)(-0.25) * (k_uy + m_uy) * (k_Bz + m_Bz)
                                  - (REAL)(-0.25) * (k_uz + m_uz) * (k_By + m_By);
  }

  inline void add_ext_num_flux_z_conservative(REAL const* uk, REAL const* um, REAL* ext_num_flux) {

      REAL m_rho = um[Field_rho];
      REAL m_E   = um[Field_E];
      REAL m_Bx  = um[Field_Bx];
      REAL m_By  = um[Field_By];
      REAL m_Bz  = um[Field_Bz];
      REAL m_ux  = um[Field_ux];
      REAL m_uy  = um[Field_uy];
      REAL m_uz  = um[Field_uz];
      REAL m_p   = um[Field_p];

      REAL m_B2 = m_Bx*m_Bx + m_By*m_By + m_Bz*m_Bz;

      REAL k_rho = uk[Field_rho];
      REAL k_E   = uk[Field_E];
      REAL k_Bx  = uk[Field_Bx];
      REAL k_By  = uk[Field_By];
      REAL k_Bz  = uk[Field_Bz];
      REAL k_ux  = uk[Field_ux];
      REAL k_uy  = uk[Field_uy];
      REAL k_uz  = uk[Field_uz];
      REAL k_p   = uk[Field_p];

      REAL k_B2 = k_Bx*k_Bx + k_By*k_By + k_Bz*k_Bz;

      ext_num_flux[Field_rho]    += (REAL)(-0.25) * (k_rho + m_rho) * (k_uz + m_uz);
      ext_num_flux[Field_rho_ux] += (REAL)(-0.25) * (k_rho*k_ux + m_rho*m_ux) * (k_uz + m_uz)
                                  - (REAL)(-0.25) * (k_Bz + m_Bz) * (k_Bx + m_Bx);
      ext_num_flux[Field_rho_uy] += (REAL)(-0.25) * (k_rho*k_uy + m_rho*m_uy) * (k_uz + m_uz)
                                  - (REAL)(-0.25) * (k_Bz + m_Bz) * (k_By + m_By);
      ext_num_flux[Field_rho_uz] += (REAL)(-0.25) * (k_rho*k_uz + m_rho*m_uz) * (k_uz + m_uz)
                                  + (REAL)(-0.5) * (k_p + (REAL)(0.5) * k_B2)
                                  - (REAL)(-0.25) * (k_Bz + m_Bz) * (k_Bz + m_Bz);
      ext_num_flux[Field_E]      += (REAL)(-0.25) * (k_E + k_p + (REAL)(0.5)*k_B2 + m_E + m_p + (REAL)(0.5)*m_B2) * (k_uz + m_uz)
                                  - (REAL)(-0.5) * (k_Bz*(k_ux*k_Bx + k_uy*k_By + k_uz*k_Bz));
      ext_num_flux[Field_Bx]     += (REAL)(-0.25) * (k_uz + m_uz) * (k_Bx + m_Bx)
                                  - (REAL)(-0.25) * (k_ux + m_ux) * (k_Bz + m_Bz);
      ext_num_flux[Field_By]     += (REAL)(-0.25) * (k_uz + m_uz) * (k_By + m_By)
                                  - (REAL)(-0.25) * (k_uy + m_uy) * (k_Bz + m_Bz);
      ext_num_flux[Field_Bz]     += (REAL)(-0.5) * (0);
  }

#else

  #error "Error in ideal_MHD.cl: No conservative discretisation specified!"

#endif


// optional nonconservative source term
#ifdef USE_SOURCE_NONE

  inline void add_ext_num_flux_x_source(REAL const* uk, REAL const* um, REAL* ext_num_flux) {
    // nothing to to
  }

  inline void add_ext_num_flux_y_source(REAL const* uk, REAL const* um, REAL* ext_num_flux) {
    // nothing to to
  }

  inline void add_ext_num_flux_z_source(REAL const* uk, REAL const* um, REAL* ext_num_flux) {
    // nothing to to
  }

#elif defined USE_SOURCE_GODUNOV_CENTRAL

  inline void add_ext_num_flux_x_source(REAL const* uk, REAL const* um, REAL* ext_num_flux) {

      REAL m_Bx  = um[Field_Bx];
      REAL m_By  = um[Field_By];
      REAL m_Bz  = um[Field_Bz];
      REAL m_ux  = um[Field_ux];
      REAL m_uy  = um[Field_uy];
      REAL m_uz  = um[Field_uz];

      REAL k_Bx  = uk[Field_Bx];
      REAL k_By  = uk[Field_By];
      REAL k_Bz  = uk[Field_Bz];

      ext_num_flux[Field_rho_ux] += (REAL)(-0.5) * (m_Bx) * (k_Bx - m_Bx);
      ext_num_flux[Field_rho_uy] += (REAL)(-0.5) * (m_By) * (k_Bx - m_Bx);
      ext_num_flux[Field_rho_uz] += (REAL)(-0.5) * (m_Bz) * (k_Bx - m_Bx);
      ext_num_flux[Field_E]      += (REAL)(-0.5) * (m_Bx*m_ux + m_By*m_uy + m_Bz*m_uz) * (k_Bx - m_Bx);
      ext_num_flux[Field_Bx]     += (REAL)(-0.5) * (m_ux) * (k_Bx - m_Bx);
      ext_num_flux[Field_By]     += (REAL)(-0.5) * (m_uy) * (k_Bx - m_Bx);
      ext_num_flux[Field_Bz]     += (REAL)(-0.5) * (m_uz) * (k_Bx - m_Bx);
  }

  inline void add_ext_num_flux_y_source(REAL const* uk, REAL const* um, REAL* ext_num_flux) {

      REAL m_Bx  = um[Field_Bx];
      REAL m_By  = um[Field_By];
      REAL m_Bz  = um[Field_Bz];
      REAL m_ux  = um[Field_ux];
      REAL m_uy  = um[Field_uy];
      REAL m_uz  = um[Field_uz];

      REAL k_Bx  = uk[Field_Bx];
      REAL k_By  = uk[Field_By];
      REAL k_Bz  = uk[Field_Bz];

      ext_num_flux[Field_rho_ux] += (REAL)(-0.5) * (m_Bx) * (k_By - m_By);
      ext_num_flux[Field_rho_uy] += (REAL)(-0.5) * (m_By) * (k_By - m_By);
      ext_num_flux[Field_rho_uz] += (REAL)(-0.5) * (m_Bz) * (k_By - m_By);
      ext_num_flux[Field_E]      += (REAL)(-0.5) * (m_Bx*m_ux + m_By*m_uy + m_Bz*m_uz) * (k_By - m_By);
      ext_num_flux[Field_Bx]     += (REAL)(-0.5) * (m_ux) * (k_By - m_By);
      ext_num_flux[Field_By]     += (REAL)(-0.5) * (m_uy) * (k_By - m_By);
      ext_num_flux[Field_Bz]     += (REAL)(-0.5) * (m_uz) * (k_By - m_By);
  }

  inline void add_ext_num_flux_z_source(REAL const* uk, REAL const* um, REAL* ext_num_flux) {

      REAL m_Bx  = um[Field_Bx];
      REAL m_By  = um[Field_By];
      REAL m_Bz  = um[Field_Bz];
      REAL m_ux  = um[Field_ux];
      REAL m_uy  = um[Field_uy];
      REAL m_uz  = um[Field_uz];

      REAL k_Bx  = uk[Field_Bx];
      REAL k_By  = uk[Field_By];
      REAL k_Bz  = uk[Field_Bz];

      ext_num_flux[Field_rho_ux] += (REAL)(-0.5) * (m_Bx) * (k_Bz - m_Bz);
      ext_num_flux[Field_rho_uy] += (REAL)(-0.5) * (m_By) * (k_Bz - m_Bz);
      ext_num_flux[Field_rho_uz] += (REAL)(-0.5) * (m_Bz) * (k_Bz - m_Bz);
      ext_num_flux[Field_E]      += (REAL)(-0.5) * (m_Bx*m_ux + m_By*m_uy + m_Bz*m_uz) * (k_Bz - m_Bz);
      ext_num_flux[Field_Bx]     += (REAL)(-0.5) * (m_ux) * (k_Bz - m_Bz);
      ext_num_flux[Field_By]     += (REAL)(-0.5) * (m_uy) * (k_Bz - m_Bz);
      ext_num_flux[Field_Bz]     += (REAL)(-0.5) * (m_uz) * (k_Bz - m_Bz);
  }

#else

  #error "Error in ideal_MHD.cl: No source discretisation specified!"

#endif



inline void compute_ext_num_flux_x(REAL const* uk, REAL const* um, REAL* ext_num_flux) {

  add_ext_num_flux_x_conservative(uk, um, ext_num_flux);
  add_ext_num_flux_x_source(uk, um, ext_num_flux);
}

inline void compute_ext_num_flux_y(REAL const* uk, REAL const* um, REAL* ext_num_flux) {

  add_ext_num_flux_y_conservative(uk, um, ext_num_flux);
  add_ext_num_flux_y_source(uk, um, ext_num_flux);
}

inline void compute_ext_num_flux_z(REAL const* uk, REAL const* um, REAL* ext_num_flux) {

  add_ext_num_flux_z_conservative(uk, um, ext_num_flux);
  add_ext_num_flux_z_source(uk, um, ext_num_flux);
}

//--------------------------------------------------------------------------------------------------
// Field initialisation
//--------------------------------------------------------------------------------------------------

inline void init_fields(uint ix, uint iy, uint iz, global REAL* u) {

  REAL rho = rho_initial(ix, iy, iz, (REAL)(0));
  REAL p = p_initial(ix, iy, iz, (REAL)(0));
  REAL4 uinit = u_initial(ix, iy, iz, (REAL)(0));
  REAL4 binit = b_initial(ix, iy, iz, (REAL)(0));

  REAL um[NUM_TOTAL_VARS] = {0};

  um[Field_rho] = rho;
  um[Field_rho_ux] = rho*uinit.x;
  um[Field_rho_uy] = rho*uinit.y;
  um[Field_rho_uz] = rho*uinit.z;
  um[Field_Bx] = binit.x;
  um[Field_By] = binit.y;
  um[Field_Bz] = binit.z;
  um[Field_E] = compute_energy(p, rho, uinit.x, uinit.y, uinit.z, binit.x, binit.y, binit.z);

  for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
    set_field_component(ix, iy, iz, i, u, um[i]);
  }
}


//--------------------------------------------------------------------------------------------------
// Surface terms
//--------------------------------------------------------------------------------------------------

#ifndef USE_PERIODIC

#define sq(x) ((x)*(x))
#define PLUS(x) (fmax(x, 0.0))

inline void calc_flux_f(const REAL *u, REAL *flux) {
  REAL BB = u[Field_Bx]*u[Field_Bx] + u[Field_By]*u[Field_By] + u[Field_Bz]*u[Field_Bz];
  flux[Field_rho] = u[Field_rho_ux];
  flux[Field_rho_ux] = u[Field_rho_ux]*u[Field_ux] + u[Field_p] + 0.5 * BB - u[Field_Bx]*u[Field_Bx];
  flux[Field_rho_uy] = u[Field_rho_ux] * u[Field_uy] - u[Field_Bx]*u[Field_By];
  flux[Field_rho_uz] = u[Field_rho_ux] * u[Field_uz] - u[Field_Bx]*u[Field_Bz];
  flux[Field_E] = u[Field_ux]*(u[Field_E] + u[Field_p] + 0.5*BB)-u[Field_Bx]*(u[Field_Bx]*u[Field_ux] + u[Field_By]*u[Field_uy] + u[Field_Bz]*u[Field_uz]);
  flux[Field_Bx] = 0;
  flux[Field_By] = u[Field_ux] * u[Field_By] - u[Field_uy]*u[Field_Bx];
  flux[Field_Bz] = u[Field_ux] * u[Field_Bz]- u[Field_uz]*u[Field_Bx];
  return;
}

inline void calc_flux_g(const REAL *u, REAL *flux){
  REAL BB = u[Field_Bx]*u[Field_Bx] + u[Field_By]*u[Field_By] + u[Field_Bz]*u[Field_Bz];
  flux[Field_rho] = u[Field_rho_uy];
  flux[Field_rho_ux] = u[Field_rho_uy] * u[Field_ux] - u[Field_By]*u[Field_Bx];
  flux[Field_rho_uy] = u[Field_rho_uy]*u[Field_uy] + u[Field_p] + 0.5 * BB - u[Field_By]*u[Field_By];
  flux[Field_rho_uz] = u[Field_rho_uy] * u[Field_uz] - u[Field_By]*u[Field_Bz];
  flux[Field_E] = u[Field_uy]*(u[Field_E] + u[Field_p] + 0.5*BB)-u[Field_By]*(u[Field_Bx]*u[Field_ux] + u[Field_By]*u[Field_uy] + u[Field_Bz]*u[Field_uz]);
  flux[Field_Bx] = u[Field_uy] * u[Field_Bx]- u[Field_ux]*u[Field_By];
  flux[Field_By] = 0;
  flux[Field_Bz] = u[Field_uy] * u[Field_Bz] - u[Field_uz]*u[Field_By];
  return;
}

inline void calc_flux_h(const REAL *u, REAL *flux){
  REAL BB = u[Field_Bx]*u[Field_Bx] + u[Field_By]*u[Field_By] + u[Field_Bz]*u[Field_Bz];
  flux[Field_rho] = u[Field_rho_uz];
  flux[Field_rho_ux] = u[Field_rho_uz] * u[Field_ux] - u[Field_Bz]*u[Field_Bx];
  flux[Field_rho_uy] = u[Field_rho_uz] * u[Field_uy] - u[Field_Bz]*u[Field_By];
  flux[Field_rho_uz] = u[Field_rho_uz] * u[Field_uz] + u[Field_p] + 0.5 * BB - u[Field_Bz]*u[Field_Bz];
  flux[Field_E] = u[Field_uz]*(u[Field_E] + u[Field_p] + 0.5*BB)-u[Field_Bz]*(u[Field_Bx]*u[Field_ux] + u[Field_By]*u[Field_uy] + u[Field_Bz]*u[Field_uz]);
  flux[Field_Bx] = u[Field_uz] * u[Field_Bx]- u[Field_ux]*u[Field_Bz];
  flux[Field_By] = u[Field_uz] * u[Field_By]- u[Field_uy]*u[Field_Bz];
  flux[Field_Bz] = 0;
  return;
}


#if defined USE_BOUNDARY_FLUX_HLL

inline void calc_num_flux(REAL al, REAL ar, REAL *ul, REAL *ur, REAL *fluxl, REAL *fluxr, REAL *num_flux){
  int i;
  if(0 < al){
    for (i = 0; i < NUM_CONSERVED_VARS; i++)
      num_flux[i] = fluxl[i];
  }
  else if(0 < ar) {
    for (i = 0; i < NUM_CONSERVED_VARS; i++)
      num_flux[i] = ((ar*fluxl[i]-al*fluxr[i]) + ar*al*(ur[i]-ul[i]))/(ar-al);
  }
  else {
    for (i = 0; i < NUM_CONSERVED_VARS; i++)
      num_flux[i] = fluxr[i];
  }
}


inline void calc_hll_speeds(REAL* ul, REAL* ur, REAL *cl, REAL* cr, int dir) {
#ifdef USE_HLL_SPEED_BKW
// Speeds from
// A multiwave approximate Riemann solver for ideal MHD based on relaxation I and II
// Bouchut, Klingenberg and Waagan (2007 and 2012)

#define alpha (0.5*(GAMMA+1)) //from page 9

  REAL rho_l = ul[Field_rho];
  REAL P_l = ul[Field_p];
  REAL rho_r = ur[Field_rho];
  REAL P_r = ur[Field_p];
  REAL u_l;
  REAL u_r;
  REAL Bx_l;
  REAL Bx_r;
  REAL Bsq_l;
  REAL Bsq_r;

  switch(dir) {
    case 0:
      u_l = ul[Field_ux];
      Bx_l = ul[Field_Bx];
      Bsq_l = ul[Field_By] * ul[Field_By] + ul[Field_Bz] * ul[Field_Bz];
      u_r = ur[Field_ux];
      Bx_r = ur[Field_Bx];
      Bsq_r = ur[Field_By] * ur[Field_By] + ur[Field_Bz] * ur[Field_Bz];
      break;
    case 1:
      u_l = ul[Field_uy];
      Bx_l = ul[Field_By];
      Bsq_l = ul[Field_Bx] * ul[Field_Bx] + ul[Field_Bz] * ul[Field_Bz];
      u_r = ur[Field_uy];
      Bx_r = ur[Field_By];
      Bsq_r = ur[Field_Bx] * ur[Field_Bx] + ur[Field_Bz] * ur[Field_Bz];
      break;
    case 2:
      u_l = ul[Field_uz];
      Bx_l = ul[Field_Bz];
      Bsq_l = ul[Field_Bx] * ul[Field_Bx] + ul[Field_By] * ul[Field_By];
      u_r = ur[Field_uz];
      Bx_r = ur[Field_Bz];
      Bsq_r = ur[Field_Bx] * ur[Field_Bx] + ur[Field_By] * ur[Field_By];
      break;
      // Default?
  }


  REAL dp_drho_cs_l = GAMMA*P_l/rho_l;
  REAL dp_drho_cs_r = GAMMA*P_r/rho_r;
  REAL pi_l = P_l + 0.5 * Bsq_l - 0.5 * Bx_l * Bx_l;
  REAL pi_r = P_r + 0.5 * Bsq_r - 0.5 * Bx_r * Bx_r;

  //fast MHD speeds (eq. 3.14)
  REAL a_ql = sqrt(0.5 * (dp_drho_cs_l + (sq(Bx_l) + Bsq_l) / rho_l + sqrt(sq(dp_drho_cs_l + (sq(Bx_l) + Bsq_l) / rho_l) - 4.0 * dp_drho_cs_l * sq(Bx_l) / rho_l)));
  REAL a_qr = sqrt(0.5 * (dp_drho_cs_r + (sq(Bx_r) + Bsq_r) / rho_r + sqrt(sq(dp_drho_cs_r + (sq(Bx_r) + Bsq_r) / rho_r) - 4.0 * dp_drho_cs_r * sq(Bx_r) / rho_r)));

  // eq. 3.17, 3.21
  REAL X_l = (PLUS(u_l - u_r) +  PLUS(pi_r - pi_l) / (rho_l * a_ql + rho_r * a_qr)) / a_ql;
  REAL x_l = 1 - X_l / (1.0 + alpha * X_l);

  // eq. 3.34
  REAL a_0l = sqrt(0.5 * (dp_drho_cs_l + (sq(Bx_l) + Bsq_l) / (x_l * rho_l) + sqrt(sq(dp_drho_cs_l + (sq(Bx_l) + Bsq_l) / (x_l * rho_l) - 4.0 * dp_drho_cs_l*sq(Bx_l) / (rho_l*x_l)))));

  // eq. 3.35
  REAL X_r = (PLUS(u_l - u_r) + PLUS(pi_l - pi_r) / (rho_l * a_ql + rho_r * a_qr)) / a_qr;
  REAL x_r = 1 - X_r / (1.0 + alpha * X_r);

  // eq 3.36
  REAL a_0r = sqrt(0.5 * (dp_drho_cs_r + (sq(Bx_r) + Bsq_r) / (x_r * rho_r) + sqrt(sq(dp_drho_cs_r + (sq(Bx_r) + Bsq_r) / (x_r * rho_r) - 4.0 * dp_drho_cs_r * sq(Bx_r) / (rho_r * x_r)))));

  // eq. 3.13
  *cl = -rho_l * a_0l - alpha * rho_l * (PLUS(u_l - u_r) + PLUS(pi_r - pi_l) / (rho_l * a_ql + rho_r * a_qr));
  *cr = rho_r * a_0r + alpha * rho_r * (PLUS(u_l - u_r) + PLUS(pi_l - pi_r) / (rho_l * a_ql + rho_r * a_qr));
  return;

#elif defined USE_HLL_SPEED_KUSANO
  // Speeds for the HLL solver are from equation (12) of
  // A multi-state HLL approximate Riemann solver for ideal magnetohydrodynamics
  // Miyoshi and Kusano (2005)

  REAL BBSQl = ul[Field_Bx] * ul[Field_Bx] + ul[Field_By] * ul[Field_By] + ul[Field_Bz] * ul[Field_Bz];
  REAL BBSQr = ur[Field_Bx] * ur[Field_Bx] + ur[Field_By] * ur[Field_By] + ur[Field_Bz] * ur[Field_Bz];

  REAL cflx, cfrx, cfly, cfry, cflz, cfrz;

  switch(dir) {
    case 0:
      // X direction
      cflx = sqrt((GAMMA * ul[Field_p] + BBSQl + sqrt(sq(GAMMA * ul[Field_p] + BBSQl) - 4 * GAMMA * ul[Field_p] * ul[Field_Bx] * ul[Field_Bx]))/(2 * ul[Field_rho]));
      cfrx = sqrt((GAMMA * ur[Field_p] + BBSQr + sqrt(sq(GAMMA * ur[Field_p] + BBSQr) - 4 * GAMMA * ur[Field_p] * ur[Field_Bx] * ur[Field_Bx]))/(2 * ur[Field_rho]));
      *cl = fmin(ul[Field_ux] - cflx, ur[Field_ux] - cfrx);
      *cr = fmax(ul[Field_ux] + cflx, ur[Field_ux] + cfrx);
      break;
    case 1:
      // Y direction
      cfly = sqrt((GAMMA * ul[Field_p] + BBSQl + sqrt(sq(GAMMA * ul[Field_p] + BBSQl) - 4 * GAMMA * ul[Field_p] * ul[Field_By] * ul[Field_By]))/(2 * ul[Field_rho]));
      cfry = sqrt((GAMMA * ur[Field_p] + BBSQr + sqrt(sq(GAMMA * ur[Field_p] + BBSQr) - 4 * GAMMA * ur[Field_p] * ur[Field_By] * ur[Field_By]))/(2 * ur[Field_rho]));
      *cl = fmin(ul[Field_uy] - cfly, ur[Field_uy] - cfry);
      *cr = fmax(ul[Field_uy] + cfly, ur[Field_uy] + cfry);
      break;
    case 2:
      // Z direction
      cflz = sqrt((GAMMA * ul[Field_p] + BBSQl + sqrt(sq(GAMMA * ul[Field_p] + BBSQl) - 4 * GAMMA * ul[Field_p] * ul[Field_Bz] * ul[Field_Bz]))/(2 * ul[Field_rho]));
      cfrz = sqrt((GAMMA * ur[Field_p] + BBSQr + sqrt(sq(GAMMA * ur[Field_p] + BBSQr) - 4 * GAMMA * ur[Field_p] * ur[Field_Bz] * ur[Field_Bz]))/(2 * ur[Field_rho]));
      *cl = fmin(ul[Field_uz] - cflz, ur[Field_uz] - cfrz);
      *cr = fmax(ul[Field_uz] + cflz, ur[Field_uz] + cfrz);
      break;
  }// Default?

#else
  #error "Define speeds to use at boundaries. Options: USE_VR_BKW / USE_VR_KUSANO."
#endif //USE_HLL_SPEED_KUSANO
}


#endif //USE_BOUNDARY_FLUX_HLL
#endif //USE_PERIODIC

inline void add_surface_terms(REAL time, uint ix, uint iy, uint iz, const global REAL *u, REAL *du_dt) {

  // For periodic boundary conditions and a single block, no surface term has to be used.
#ifndef USE_PERIODIC
#ifdef USE_BOUNDARY_FLUX_HLL
  int i;
  REAL um[NUM_TOTAL_VARS] = {0.0};
  get_field(ix, iy, iz, 0, 0, 0, u, um);
  REAL4 b_bound = b_boundary(ix, iy, iz, time);
  REAL4 u_bound = u_boundary(ix, iy, iz, time);
  REAL rho_bound = rho_boundary(ix, iy, iz, time); //(REAL)1;
  REAL p_bound = p_boundary(ix, iy, iz, time);
  REAL E_bound = compute_energy(p_bound, rho_bound, u_bound.x, u_bound.y, u_bound.z, b_bound.x, b_bound.y, b_bound.z);


  REAL ub[NUM_TOTAL_VARS] = {0.0};
  ub[Field_ux] = u_bound.x;
  ub[Field_uy] = u_bound.y;
  ub[Field_uz] = u_bound.z;
  ub[Field_rho_ux] = u_bound.x * rho_bound;
  ub[Field_rho_uy] = u_bound.y * rho_bound;
  ub[Field_rho_uz] = u_bound.z * rho_bound;
  ub[Field_rho] = rho_bound;
  ub[Field_p] = p_bound;
  ub[Field_E] = E_bound;
  ub[Field_Bx] = b_bound.x;
  ub[Field_By] = b_bound.y;
  ub[Field_Bz] = b_bound.z;

  REAL flux[NUM_CONSERVED_VARS] = {0.0};
  REAL fluxm[NUM_CONSERVED_VARS] = {0.0};
  REAL fluxb[NUM_CONSERVED_VARS] = {0.0};

  REAL alx, arx, aly, ary, alz, arz;

  if (check_bound_xr(ix, 1)) {
    calc_hll_speeds(um, ub, &alx, &arx, 0);
    calc_flux_f(um, fluxm);
    calc_flux_f(ub, fluxb);
    calc_num_flux(alx, arx, um, ub, fluxm, fluxb, flux);
    for(i = 0; i < NUM_CONSERVED_VARS; i++)
      du_dt[i] -= (REAL)(M_INV[0] / DX) * (flux[i] - fluxm[i]);
  }
  else if(check_bound_l(ix,1)) {
    calc_hll_speeds(ub, um, &alx, &arx, 0);
    calc_flux_f(um, fluxm);
    calc_flux_f(ub, fluxb);
    calc_num_flux(alx, arx, ub, um, fluxb, fluxm, flux);
    for(i = 0; i < NUM_CONSERVED_VARS; i++)
      du_dt[i] += (REAL)(M_INV[0] / DX) * (flux[i] - fluxm[i]);
  }

  if (check_bound_yr(iy, 1)) {
    calc_hll_speeds(um, ub, &aly, &ary, 1);
    calc_flux_g(um, fluxm);
    calc_flux_g(ub, fluxb);
    calc_num_flux(aly, ary, um, ub, fluxm, fluxb, flux);
    for(i = 0; i < NUM_CONSERVED_VARS; i++)
      du_dt[i] -= (REAL)(M_INV[0] / DY) * (flux[i] - fluxm[i]);
  }
  else if(check_bound_l(iy,1)) {
    calc_hll_speeds(ub, um, &aly, &ary, 1);
    calc_flux_g(um, fluxm);
    calc_flux_g(ub, fluxb);
    calc_num_flux(aly, ary, ub, um, fluxb, fluxm, flux);
    for(i = 0; i < NUM_CONSERVED_VARS; i++)
      du_dt[i] += (REAL)(M_INV[0] / DY) * (flux[i] - fluxm[i]);
  }

  if (check_bound_zr(iz, 1)) {
    calc_hll_speeds(um, ub, &alz, &arz, 2);
    calc_flux_h(um, fluxm);
    calc_flux_h(ub, fluxb);
    calc_num_flux(alz, arz, um, ub, fluxm, fluxb, flux);
    for(i = 0; i < NUM_CONSERVED_VARS; i++)
      du_dt[i] -= (REAL)(M_INV[0] / DZ) * (flux[i] - fluxm[i]);
  }
  else if(check_bound_l(iz,1)) {
    calc_hll_speeds(ub, um, &alz, &arz, 2);
    calc_flux_h(um, fluxm);
    calc_flux_h(ub, fluxb);
    calc_num_flux(alz,arz, ub, um, fluxb, fluxm, flux);
    for(i = 0; i < NUM_CONSERVED_VARS; i++)
      du_dt[i] += (REAL)(M_INV[0] / DZ) * (flux[i] - fluxm[i]);
  }

#else
  #error "Error in include_physics/ideal_MHD.h:Non-periodic boundarys are used but no boundary flux is specified!"

#endif // USE_BOUNDARY_FLUX_HLL
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
  REAL Bx     = um[Field_Bx];
  REAL By     = um[Field_By];
  REAL Bz     = um[Field_Bz];

  // velocity
  REAL ux = rho_ux / rho;
  REAL uy = rho_uy / rho;
  REAL uz = rho_uz / rho;
  set_field_component(ix, iy, iz, Field_ux, u, ux);
  set_field_component(ix, iy, iz, Field_uy, u, uy);
  set_field_component(ix, iy, iz, Field_uz, u, uz);

  // pressure
  REAL p = compute_pressure(rho, ux, uy, uz, E, Bx, By, Bz);
  set_field_component(ix, iy, iz, Field_p, u, p);
}

//--------------------------------------------------------------------------------------------------
// Computation of analytical solution for all conserved variables
//---------------------------------------------------------------------------------------------------

inline void analytical_solution(uint ix, uint iy, uint iz, global REAL *u, REAL time) {

  REAL rho = (REAL)(1);
  REAL p = (REAL)(1);
  REAL4 B =  b_analytical(ix, iy, iz, time);
  REAL4 v =  u_analytical(ix, iy, iz, time);
  REAL E = compute_energy(p, rho, v.x, v.y, v.z, B.x, B.y, B.z);

  set_field_component(ix, iy, iz, Field_rho, u, rho);

  set_field_component(ix, iy, iz, Field_rho_ux, u, rho*v.x);
  set_field_component(ix, iy, iz, Field_rho_uy, u, rho*v.y);
  set_field_component(ix, iy, iz, Field_rho_uz, u, rho*v.z);

  set_field_component(ix, iy, iz, Field_E, u, E);

  set_field_component(ix, iy, iz, Field_Bx, u, B.x);
  set_field_component(ix, iy, iz, Field_By, u, B.y);
  set_field_component(ix, iy, iz, Field_Bz, u, B.z);
}

#endif // IDEAL_MHD
