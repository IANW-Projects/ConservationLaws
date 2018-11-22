// This project is licensed under the terms of the Creative Commons CC BY-NC-ND 4.0 license.

#ifndef IDEAL_MHD
#define IDEAL_MHD


// ideal gas constant
REAL constant GAMMA = (REAL)(1.4);

// paramters of the Alfven wave test case
REAL constant CONST_A = 0.1;
REAL constant CONST_alpha = 0.5235987755982988; // pi/6


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
// Functions for field initialisation
//--------------------------------------------------------------------------------------------------

/*
Alfven wave testcase used in section 8.6 of
Skew-Symmetric Splitting for Multiscale Gas Dynamics and MHD Turbulence Flows.
Björn Sjögreen, Helen C. Yee, Dmitry Kotov, 2017.
*/

/*
Analytical solution of the velocity field.
*/
inline REAL4 u_analytical(uint ix, uint iy, uint iz, REAL time) {

	REAL x = (REAL)XMIN + ix*(REAL)DX;
	REAL y = (REAL)YMIN + iy*(REAL)DY;
	REAL z = (REAL)ZMIN + iz*(REAL)DZ;

  REAL cos_alpha = cos(CONST_alpha);
  REAL sin_alpha = sin(CONST_alpha);

	REAL phi = (REAL)(6.283185307179586) * (x*cos_alpha + y*sin_alpha + time);
  REAL cos_phi = cos(phi);
  REAL sin_phi = sin(phi);

  REAL u1 = - CONST_A * sin_alpha * sin_phi;
  REAL u2 = CONST_A * cos_alpha * sin_phi;
  REAL u3 = CONST_A * cos_phi;

	return (REAL4) {u1, u2, u3, (REAL)(0)};
}

/*
Analytical solution of the magnetic field.
*/
inline REAL4 b_analytical(uint ix, uint iy, uint iz, REAL time) {

	REAL x = (REAL)XMIN + ix*(REAL)DX;
	REAL y = (REAL)YMIN + iy*(REAL)DY;
	REAL z = (REAL)ZMIN + iz*(REAL)DZ;

  REAL cos_alpha = cos(CONST_alpha);
  REAL sin_alpha = sin(CONST_alpha);

	REAL phi = (REAL)(6.283185307179586) * (x*cos_alpha + y*sin_alpha + time);
  REAL cos_phi = cos(phi);
  REAL sin_phi = sin(phi);

  REAL B1 = cos_alpha - CONST_A * sin_alpha * sin_phi;
  REAL B2 = sin_alpha + CONST_A * cos_alpha * sin_phi;
  REAL B3 = CONST_A * cos_phi;

	return (REAL4) {B1, B2, B3, (REAL)(0)};
}


inline void init_fields(uint ix, uint iy, uint iz, global REAL* u) {

  REAL rho = (REAL)(1);
  REAL p = (REAL)(1);
  REAL4 uinit = u_analytical(ix, iy, iz, (REAL)(0));
  REAL4 binit = b_analytical(ix, iy, iz, (REAL)(0));

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


inline void add_surface_terms(REAL time, uint ix, uint iy, uint iz, global REAL *u, REAL *du_dt) {

  // For periodic boundary conditions and a single block, no surface term has to be used.
#ifndef USE_PERIODIC
  #error "Error in ideal_MHD.cl: Surface terms not implemented for nonperiodic boundaries!"
#endif // USE_PERIODIC*/

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



#endif // IDEAL_MHD
