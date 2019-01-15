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


inline void compute_flux_x(REAL* u, REAL* flux) {
  flux[Field_rho]    = u[Field_rho_ux];
  flux[Field_rho_ux] = u[Field_rho_ux] * u[Field_ux] + u[Field_p];
  flux[Field_rho_uy] = u[Field_rho_ux] * u[Field_uy];
  flux[Field_rho_uz] = u[Field_rho_ux] * u[Field_uz];
  flux[Field_E]      = (u[Field_E] + u[Field_p]) * u[Field_ux];
}

inline void compute_flux_y(REAL* u, REAL* flux) {
  flux[Field_rho]    = u[Field_rho_uy];
  flux[Field_rho_ux] = u[Field_rho_uy] * u[Field_ux];
  flux[Field_rho_uy] = u[Field_rho_uy] * u[Field_uy] + u[Field_p];
  flux[Field_rho_uz] = u[Field_rho_uy] * u[Field_uz];
  flux[Field_E]      = (u[Field_E] + u[Field_p]) * u[Field_uy];
}

inline void compute_flux_z(REAL* u, REAL* flux) {
  flux[Field_rho]    = u[Field_rho_uz];
  flux[Field_rho_ux] = u[Field_rho_uz] * u[Field_ux];
  flux[Field_rho_uy] = u[Field_rho_uz] * u[Field_uy];
  flux[Field_rho_uz] = u[Field_rho_uz] * u[Field_uz] + u[Field_p];
  flux[Field_E]      = (u[Field_E] + u[Field_p]) * u[Field_uz];
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
/* Numerical fluxes of
@article{ducros2000high,
  title={High-Order Fluxes for Conservative Skew-Symmetric-Like Schemes in
         Structured Meshes: {A}pplication to Compressible Flows},
  author={Ducros, F and Laporte, F and Souleres, T and Guinot, V and Moinat, P
          and Caruelle, B},
  journal={Journal of Computational Physics},
  volume={161},
  number={1},
  pages={114--139},
  year={2000},
  publisher={Elsevier},
  doi={10.1006/jcph.2000.6492}
}
*/

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

#elif defined USE_FLUX_KennedyGruber
/* Kinetic energy preserving numerical fluxes of
@article{kennedy2008reduced,
  title={Reduced aliasing formulations of the convective terms within the
        {N}avier--{S}tokes equations for a compressible fluid},
  author={Kennedy, Christopher A and Gruber, Andrea},
  journal={Journal of Computational Physics},
  volume={227},
  number={3},
  pages={1676--1700},
  year={2008},
  publisher={Elsevier},
  doi={10.1016/j.jcp.2007.09.020}
}
*/

  inline void compute_ext_num_flux_x(REAL const* uk, REAL const* um, REAL* ext_num_flux) {

      REAL m_rho = um[Field_rho];
      REAL m_E   = um[Field_E];
      REAL m_ux  = um[Field_ux];
      REAL m_uy  = um[Field_uy];
      REAL m_uz  = um[Field_uz];
      REAL m_p   = um[Field_p];
      REAL m_e   = m_E / m_rho;

      REAL k_rho = uk[Field_rho];
      REAL k_E   = uk[Field_E];
      REAL k_ux  = uk[Field_ux];
      REAL k_uy  = uk[Field_uy];
      REAL k_uz  = uk[Field_uz];
      REAL k_p   = uk[Field_p];
      REAL k_e   = k_E / k_rho;

      ext_num_flux[Field_rho]    = (REAL)(-0.25) * (k_rho + m_rho) * (k_ux + m_ux);
      ext_num_flux[Field_rho_ux] = (REAL)(-0.125) * (k_rho + m_rho) * (k_ux + m_ux) * (k_ux + m_ux) + (REAL)(-0.5) * (k_p + m_p);
      ext_num_flux[Field_rho_uy] = (REAL)(-0.125) * (k_rho + m_rho) * (k_ux + m_ux) * (k_uy + m_uy);
      ext_num_flux[Field_rho_uz] = (REAL)(-0.125) * (k_rho + m_rho) * (k_ux + m_ux) * (k_uz + m_uz);
      ext_num_flux[Field_E]      = (REAL)(-0.25) * ((REAL)(0.5) * (k_rho + m_rho) * (k_e + m_e) + (k_p + m_p)) * (k_ux + m_ux);
  }

  inline void compute_ext_num_flux_y(REAL const* uk, REAL const* um, REAL* ext_num_flux) {

      REAL m_rho = um[Field_rho];
      REAL m_E   = um[Field_E];
      REAL m_ux  = um[Field_ux];
      REAL m_uy  = um[Field_uy];
      REAL m_uz  = um[Field_uz];
      REAL m_p   = um[Field_p];
      REAL m_e   = m_E / m_rho;

      REAL k_rho = uk[Field_rho];
      REAL k_E   = uk[Field_E];
      REAL k_ux  = uk[Field_ux];
      REAL k_uy  = uk[Field_uy];
      REAL k_uz  = uk[Field_uz];
      REAL k_p   = uk[Field_p];
      REAL k_e   = k_E / k_rho;

      ext_num_flux[Field_rho]    = (REAL)(-0.25) * (k_rho + m_rho) * (k_uy + m_uy);
      ext_num_flux[Field_rho_ux] = (REAL)(-0.125) * (k_rho + m_rho) * (k_uy + m_uy) * (k_ux + m_ux);
      ext_num_flux[Field_rho_uy] = (REAL)(-0.125) * (k_rho + m_rho) * (k_uy + m_uy) * (k_uy + m_uy) + (REAL)(-0.5) * (k_p + m_p);
      ext_num_flux[Field_rho_uz] = (REAL)(-0.125) * (k_rho + m_rho) * (k_uy + m_uy) * (k_uz + m_uz);
      ext_num_flux[Field_E]      = (REAL)(-0.25) * ((REAL)(0.5) * (k_rho + m_rho) * (k_e + m_e) + (k_p + m_p)) * (k_uy + m_uy);
  }

  inline void compute_ext_num_flux_z(REAL const* uk, REAL const* um, REAL* ext_num_flux) {

      REAL m_rho = um[Field_rho];
      REAL m_E   = um[Field_E];
      REAL m_ux  = um[Field_ux];
      REAL m_uy  = um[Field_uy];
      REAL m_uz  = um[Field_uz];
      REAL m_p   = um[Field_p];
      REAL m_e   = m_E / m_rho;

      REAL k_rho = uk[Field_rho];
      REAL k_E   = uk[Field_E];
      REAL k_ux  = uk[Field_ux];
      REAL k_uy  = uk[Field_uy];
      REAL k_uz  = uk[Field_uz];
      REAL k_p   = uk[Field_p];
      REAL k_e   = k_E / k_rho;

      ext_num_flux[Field_rho]    = (REAL)(-0.25) * (k_rho + m_rho) * (k_uz + m_uz);
      ext_num_flux[Field_rho_ux] = (REAL)(-0.125) * (k_rho + m_rho) * (k_uz + m_uz) * (k_ux + m_ux);
      ext_num_flux[Field_rho_uy] = (REAL)(-0.125) * (k_rho + m_rho) * (k_uz + m_uz) * (k_uy + m_uy);
      ext_num_flux[Field_rho_uz] = (REAL)(-0.125) * (k_rho + m_rho) * (k_uz + m_uz) * (k_uz + m_uz) + (REAL)(-0.5) * (k_p + m_p);
      ext_num_flux[Field_E]      = (REAL)(-0.25) * ((REAL)(0.5) * (k_rho + m_rho) * (k_e + m_e) + (k_p + m_p)) * (k_uz + m_uz);
  }

#elif defined USE_FLUX_Morinishi
/* Kinetic energy preserving numerical fluxes of
@article{morinishi2010skew,
  title={Skew-symmetric form of convective terms and fully conservative finite
        difference schemes for variable density low-{M}ach number flows},
  author={Morinishi, Yohei},
  journal={Journal of Computational Physics},
  volume={229},
  number={2},
  pages={276--300},
  year={2010},
  publisher={Elsevier},
  doi={10.1016/j.jcp.2009.09.021}
}
*/

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

      ext_num_flux[Field_rho]    = (REAL)(-0.5) * (k_rho*k_ux + m_rho*m_ux);
      ext_num_flux[Field_rho_ux] = (REAL)(-0.25) * (k_rho*k_ux + m_rho*m_ux) * (k_ux + m_ux) + (REAL)(-0.5) * (k_p + m_p);
      ext_num_flux[Field_rho_uy] = (REAL)(-0.25) * (k_rho*k_ux + m_rho*m_ux) * (k_uy + m_uy);
      ext_num_flux[Field_rho_uz] = (REAL)(-0.25) * (k_rho*k_ux + m_rho*m_ux) * (k_uz + m_uz);
      ext_num_flux[Field_E]      = (REAL)(GAMMA/(GAMMA-1)) * (REAL)(-0.5) * (k_p*k_ux + m_p*m_ux)
                                  + (REAL)(-0.25) * (k_rho*k_ux*k_ux + m_rho*m_ux*m_ux) * (k_ux + m_ux)
                                  + (REAL)(-0.25) * (k_rho*k_ux*k_uy + m_rho*m_ux*m_uy) * (k_uy + m_uy)
                                  + (REAL)(-0.25) * (k_rho*k_ux*k_uz + m_rho*m_ux*m_uz) * (k_uz + m_uz)
                                  - (REAL)(-0.25) * (k_rho*k_ux*(k_ux*k_ux+k_uy*k_uy+k_uz*k_uz) + m_rho*m_ux*(m_ux*m_ux+m_uy*m_uy+m_uz*m_uz));
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

      ext_num_flux[Field_rho]    = (REAL)(-0.5) * (k_rho*k_uy + m_rho*m_uy);
      ext_num_flux[Field_rho_ux] = (REAL)(-0.25) * (k_rho*k_uy + m_rho*m_uy) * (k_ux + m_ux);
      ext_num_flux[Field_rho_uy] = (REAL)(-0.25) * (k_rho*k_uy + m_rho*m_uy) * (k_uy + m_uy) + (REAL)(-0.5) * (k_p + m_p);
      ext_num_flux[Field_rho_uz] = (REAL)(-0.25) * (k_rho*k_uy + m_rho*m_uy) * (k_uz + m_uz);
      ext_num_flux[Field_E]      = (REAL)(GAMMA/(GAMMA-1)) * (REAL)(-0.5) * (k_p*k_uy + m_p*m_uy)
                                  + (REAL)(-0.25) * (k_rho*k_uy*k_ux + m_rho*m_uy*m_ux) * (k_ux + m_ux)
                                  + (REAL)(-0.25) * (k_rho*k_uy*k_uy + m_rho*m_uy*m_uy) * (k_uy + m_uy)
                                  + (REAL)(-0.25) * (k_rho*k_uy*k_uz + m_rho*m_uy*m_uz) * (k_uz + m_uz)
                                  - (REAL)(-0.25) * (k_rho*k_uy*(k_ux*k_ux+k_uy*k_uy+k_uz*k_uz) + m_rho*m_uy*(m_ux*m_ux+m_uy*m_uy+m_uz*m_uz));
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

      ext_num_flux[Field_rho]    = (REAL)(-0.5) * (k_rho*k_uz + m_rho*m_uz);
      ext_num_flux[Field_rho_ux] = (REAL)(-0.25) * (k_rho*k_uz + m_rho*m_uz) * (k_ux + m_ux);
      ext_num_flux[Field_rho_uy] = (REAL)(-0.25) * (k_rho*k_uz + m_rho*m_uz) * (k_uy + m_uy);
      ext_num_flux[Field_rho_uz] = (REAL)(-0.25) * (k_rho*k_uz + m_rho*m_uz) * (k_uz + m_uz) + (REAL)(-0.5) * (k_p + m_p);
      ext_num_flux[Field_E]      = (REAL)(GAMMA/(GAMMA-1)) * (REAL)(-0.5) * (k_p*k_uz + m_p*m_uz)
                                  + (REAL)(-0.25) * (k_rho*k_uz*k_ux + m_rho*m_uz*m_ux) * (k_ux + m_ux)
                                  + (REAL)(-0.25) * (k_rho*k_uz*k_uy + m_rho*m_uz*m_uy) * (k_uy + m_uy)
                                  + (REAL)(-0.25) * (k_rho*k_uz*k_uz + m_rho*m_uz*m_uz) * (k_uz + m_uz)
                                  - (REAL)(-0.25) * (k_rho*k_uz*(k_ux*k_ux+k_uy*k_uy+k_uz*k_uz) + m_rho*m_uz*(m_ux*m_ux+m_uy*m_uy+m_uz*m_uz));
  }

#elif defined USE_FLUX_Pirozzoli
/* Kinetic energy preserving numerical fluxes of
@article{pirozzoli2011numerical,
  title={Numerical Methods for High-Speed Flows},
  author={Pirozzoli, Sergio},
  journal={Annual Review of Fluid Mechanics},
  volume={43},
  pages={163--194},
  year={2011},
  publisher={Annual Reviews},
  doi={10.1146/annurev-fluid-122109-160718}
}
*/

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
      ext_num_flux[Field_rho_ux] = (REAL)(-0.125) * (k_rho + m_rho) * (k_ux + m_ux) * (k_ux + m_ux) + (REAL)(-0.5) * (k_p + m_p);
      ext_num_flux[Field_rho_uy] = (REAL)(-0.125) * (k_rho + m_rho) * (k_ux + m_ux) * (k_uy + m_uy);
      ext_num_flux[Field_rho_uz] = (REAL)(-0.125) * (k_rho + m_rho) * (k_ux + m_ux) * (k_uz + m_uz);
      ext_num_flux[Field_E]      = (REAL)(-0.125) * (k_rho + m_rho) * ((k_E+k_p)/k_rho + (m_E+m_p)/m_rho) * (k_ux + m_ux);
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
      ext_num_flux[Field_rho_ux] = (REAL)(-0.125) * (k_rho + m_rho) * (k_uy + m_uy) * (k_ux + m_ux);
      ext_num_flux[Field_rho_uy] = (REAL)(-0.125) * (k_rho + m_rho) * (k_uy + m_uy) * (k_uy + m_uy) + (REAL)(-0.5) * (k_p + m_p);
      ext_num_flux[Field_rho_uz] = (REAL)(-0.125) * (k_rho + m_rho) * (k_uy + m_uy) * (k_uz + m_uz);
      ext_num_flux[Field_E]      = (REAL)(-0.125) * (k_rho + m_rho) * ((k_E+k_p)/k_rho + (m_E+m_p)/m_rho) * (k_uy + m_uy);
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
      ext_num_flux[Field_rho_ux] = (REAL)(-0.125) * (k_rho + m_rho) * (k_uz + m_uz) * (k_ux + m_ux);
      ext_num_flux[Field_rho_uy] = (REAL)(-0.125) * (k_rho + m_rho) * (k_uz + m_uz) * (k_uy + m_uy);
      ext_num_flux[Field_rho_uz] = (REAL)(-0.125) * (k_rho + m_rho) * (k_uz + m_uz) * (k_uz + m_uz) + (REAL)(-0.5) * (k_p + m_p);
      ext_num_flux[Field_E]      = (REAL)(-0.125) * (k_rho + m_rho) * ((k_E+k_p)/k_rho + (m_E+m_p)/m_rho) * (k_uz + m_uz);
  }

#elif defined USE_FLUX_KuyaTotaniKawai_QCC
/* Kinetic energy preserving numerical fluxes of
@article{kuya2018kinetic,
  title={Kinetic energy and entropy preserving schemes for compressible flows
         by split convective forms},
  author={Kuya, Yuichi and Totani, Kosuke and Kawai, Soshi},
  journal={Journal of Computational Physics},
  volume={375},
  pages={823--853},
  year={2018},
  publisher={Elsevier},
  doi={10.1016/j.jcp.2018.08.058}
}
*/

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
      ext_num_flux[Field_rho_ux] = (REAL)(-0.125) * (k_rho + m_rho) * (k_ux + m_ux) * (k_ux + m_ux) + (REAL)(-0.5) * (k_p + m_p);
      ext_num_flux[Field_rho_uy] = (REAL)(-0.125) * (k_rho + m_rho) * (k_uy + m_uy) * (k_ux + m_ux);
      ext_num_flux[Field_rho_uz] = (REAL)(-0.125) * (k_rho + m_rho) * (k_uz + m_uz) * (k_ux + m_ux);
      ext_num_flux[Field_E]      = (REAL)(-0.125) * (k_rho + m_rho) * (k_ux + m_ux)
                                                  * (k_ux*m_ux + k_uy*m_uy + k_uz*m_uz + (k_p/k_rho + m_p/m_rho)/(GAMMA-1))
                                 + (REAL)(-0.5) * (k_ux*m_p + m_ux*k_p);
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
      ext_num_flux[Field_rho_ux] = (REAL)(-0.125) * (k_rho + m_rho) * (k_ux + m_ux) * (k_uy + m_uy);
      ext_num_flux[Field_rho_uy] = (REAL)(-0.125) * (k_rho + m_rho) * (k_uy + m_uy) * (k_uy + m_uy) + (REAL)(-0.5) * (k_p + m_p);
      ext_num_flux[Field_rho_uz] = (REAL)(-0.125) * (k_rho + m_rho) * (k_uz + m_uz) * (k_uy + m_uy);
      ext_num_flux[Field_E]      = (REAL)(-0.125) * (k_rho + m_rho) * (k_uy + m_uy)
                                                  * (k_ux*m_ux + k_uy*m_uy + k_uz*m_uz + (k_p/k_rho + m_p/m_rho)/(GAMMA-1))
                                 + (REAL)(-0.5) * (k_uy*m_p + m_uy*k_p);
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
      ext_num_flux[Field_rho_ux] = (REAL)(-0.125) * (k_rho + m_rho) * (k_ux + m_ux) * (k_uz + m_uz);
      ext_num_flux[Field_rho_uy] = (REAL)(-0.125) * (k_rho + m_rho) * (k_uy + m_uy) * (k_uz + m_uz);
      ext_num_flux[Field_rho_uz] = (REAL)(-0.125) * (k_rho + m_rho) * (k_uz + m_uz) * (k_uz + m_uz) + (REAL)(-0.5) * (k_p + m_p);
      ext_num_flux[Field_E]      = (REAL)(-0.125) * (k_rho + m_rho) * (k_uz + m_uz)
                                                  * (k_ux*m_ux + k_uy*m_uy + k_uz*m_uz + (k_p/k_rho + m_p/m_rho)/(GAMMA-1))
                                 + (REAL)(-0.5) * (k_uz*m_p + m_uz*k_p);
  }

#elif defined USE_FLUX_Chandrashekar
/* Entropy conservative (and kinetic energy preserving) numerical fluxes of
@article{chandrashekar2013kinetic,
  title={Kinetic Energy Preserving and Entropy Stable Finite Volume Schemes for
        Compressible {E}uler and {N}avier-{S}tokes Equations},
  author={Chandrashekar, Praveen},
  journal={Communications in Computational Physics},
  volume={14},
  number={5},
  pages={1252--1286},
  year={2013},
  doi={10.4208/cicp.170712.010313a}
}
*/

  inline void compute_ext_num_flux_x(REAL const* uk, REAL const* um, REAL* ext_num_flux) {

      REAL m_rho = um[Field_rho];
      REAL m_E   = um[Field_E];
      REAL m_ux  = um[Field_ux];
      REAL m_uy  = um[Field_uy];
      REAL m_uz  = um[Field_uz];
      REAL m_p   = um[Field_p];
      REAL m_beta= m_rho / (2 * m_p);

      REAL k_rho = uk[Field_rho];
      REAL k_E   = uk[Field_E];
      REAL k_ux  = uk[Field_ux];
      REAL k_uy  = uk[Field_uy];
      REAL k_uz  = uk[Field_uz];
      REAL k_p   = uk[Field_p];
      REAL k_beta= k_rho / (2 * k_p);

      REAL rho     = (REAL)(0.5) * (m_rho + k_rho);
      REAL rho_log = logmean(m_rho, k_rho);
      REAL ux      = (REAL)(0.5) * (m_ux + k_ux);
      REAL uy      = (REAL)(0.5) * (m_uy + k_uy);
      REAL uz      = (REAL)(0.5) * (m_uz + k_uz);
      REAL u2      = (REAL)(0.5) * (m_ux*m_ux+m_uy*m_uy+m_uz*m_uz + k_ux*k_ux+k_uy*k_uy+k_uz*k_uz);
      REAL beta    = (REAL)(0.5) * (m_beta + k_beta);
      REAL beta_log= logmean(m_beta, k_beta);

      REAL f_rho   = rho_log*ux;
      REAL f_rhoux = ux*f_rho + rho/(2*beta);
      REAL f_rhouy = uy*f_rho;
      REAL f_rhouz = uz*f_rho;
      REAL f_E     = (REAL)(1/(2*GAMMA-2))*f_rho/beta_log - (REAL)(0.5)*u2*f_rho + ux*f_rhoux + uy*f_rhouy + uz*f_rhouz;

      ext_num_flux[Field_rho]    = - f_rho;
      ext_num_flux[Field_rho_ux] = - f_rhoux;
      ext_num_flux[Field_rho_uy] = - f_rhouy;
      ext_num_flux[Field_rho_uz] = - f_rhouz;
      ext_num_flux[Field_E]      = - f_E;
  }

  inline void compute_ext_num_flux_y(REAL const* uk, REAL const* um, REAL* ext_num_flux) {

      REAL m_rho = um[Field_rho];
      REAL m_E   = um[Field_E];
      REAL m_ux  = um[Field_ux];
      REAL m_uy  = um[Field_uy];
      REAL m_uz  = um[Field_uz];
      REAL m_p   = um[Field_p];
      REAL m_beta= m_rho / (2 * m_p);

      REAL k_rho = uk[Field_rho];
      REAL k_E   = uk[Field_E];
      REAL k_ux  = uk[Field_ux];
      REAL k_uy  = uk[Field_uy];
      REAL k_uz  = uk[Field_uz];
      REAL k_p   = uk[Field_p];
      REAL k_beta= k_rho / (2 * k_p);

      REAL rho     = (REAL)(0.5) * (m_rho + k_rho);
      REAL rho_log = logmean(m_rho, k_rho);
      REAL ux      = (REAL)(0.5) * (m_ux + k_ux);
      REAL uy      = (REAL)(0.5) * (m_uy + k_uy);
      REAL uz      = (REAL)(0.5) * (m_uz + k_uz);
      REAL u2      = (REAL)(0.5) * (m_ux*m_ux+m_uy*m_uy+m_uz*m_uz + k_ux*k_ux+k_uy*k_uy+k_uz*k_uz);
      REAL beta    = (REAL)(0.5) * (m_beta + k_beta);
      REAL beta_log= logmean(m_beta, k_beta);

      REAL f_rho   = rho_log*uy;
      REAL f_rhoux = ux*f_rho;
      REAL f_rhouy = uy*f_rho + rho/(2*beta);
      REAL f_rhouz = uz*f_rho;
      REAL f_E     = (REAL)(1/(2*GAMMA-2))*f_rho/beta_log - (REAL)(0.5)*u2*f_rho + ux*f_rhoux + uy*f_rhouy + uz*f_rhouz;

      ext_num_flux[Field_rho]    = - f_rho;
      ext_num_flux[Field_rho_ux] = - f_rhoux;
      ext_num_flux[Field_rho_uy] = - f_rhouy;
      ext_num_flux[Field_rho_uz] = - f_rhouz;
      ext_num_flux[Field_E]      = - f_E;
  }

  inline void compute_ext_num_flux_z(REAL const* uk, REAL const* um, REAL* ext_num_flux) {

      REAL m_rho = um[Field_rho];
      REAL m_E   = um[Field_E];
      REAL m_ux  = um[Field_ux];
      REAL m_uy  = um[Field_uy];
      REAL m_uz  = um[Field_uz];
      REAL m_p   = um[Field_p];
      REAL m_beta= m_rho / (2 * m_p);

      REAL k_rho = uk[Field_rho];
      REAL k_E   = uk[Field_E];
      REAL k_ux  = uk[Field_ux];
      REAL k_uy  = uk[Field_uy];
      REAL k_uz  = uk[Field_uz];
      REAL k_p   = uk[Field_p];
      REAL k_beta= k_rho / (2 * k_p);

      REAL rho     = (REAL)(0.5) * (m_rho + k_rho);
      REAL rho_log = logmean(m_rho, k_rho);
      REAL ux      = (REAL)(0.5) * (m_ux + k_ux);
      REAL uy      = (REAL)(0.5) * (m_uy + k_uy);
      REAL uz      = (REAL)(0.5) * (m_uz + k_uz);
      REAL u2      = (REAL)(0.5) * (m_ux*m_ux+m_uy*m_uy+m_uz*m_uz + k_ux*k_ux+k_uy*k_uy+k_uz*k_uz);
      REAL beta    = (REAL)(0.5) * (m_beta + k_beta);
      REAL beta_log= logmean(m_beta, k_beta);

      REAL f_rho   = rho_log*uz;
      REAL f_rhoux = ux*f_rho;
      REAL f_rhouy = uy*f_rho;
      REAL f_rhouz = uz*f_rho + rho/(2*beta);
      REAL f_E     = (REAL)(1/(2*GAMMA-2))*f_rho/beta_log - (REAL)(0.5)*u2*f_rho + ux*f_rhoux + uy*f_rhouy + uz*f_rhouz;

      ext_num_flux[Field_rho]    = - f_rho;
      ext_num_flux[Field_rho_ux] = - f_rhoux;
      ext_num_flux[Field_rho_uy] = - f_rhouy;
      ext_num_flux[Field_rho_uz] = - f_rhouz;
      ext_num_flux[Field_E]      = - f_E;
  }

#elif defined USE_FLUX_IsmailRoe
/* Entropy conservative numerical fluxes of
@article{ismail2009affordable,
  title={Affordable, entropy-consistent {E}uler flux functions {II}:
        {E}ntropy production at shocks},
  author={Ismail, Farzad and Roe, Philip L},
  journal={Journal of Computational Physics},
  volume={228},
  number={15},
  pages={5410--5436},
  year={2009},
  publisher={Elsevier},
  doi={10.1016/j.jcp.2009.04.021}
}
*/

  inline void compute_ext_num_flux_x(REAL const* uk, REAL const* um, REAL* ext_num_flux) {

      REAL m_rho = um[Field_rho];
      REAL m_E   = um[Field_E];
      REAL m_ux  = um[Field_ux];
      REAL m_uy  = um[Field_uy];
      REAL m_uz  = um[Field_uz];
      REAL m_p   = um[Field_p];
      REAL m_z1  = sqrt(m_rho / m_p);
      REAL m_z2  = m_z1 * m_ux;
      REAL m_z3  = m_z1 * m_uy;
      REAL m_z4  = m_z1 * m_uz;
      REAL m_z5  = m_z1 * m_p;

      REAL k_rho = uk[Field_rho];
      REAL k_E   = uk[Field_E];
      REAL k_ux  = uk[Field_ux];
      REAL k_uy  = uk[Field_uy];
      REAL k_uz  = uk[Field_uz];
      REAL k_p   = uk[Field_p];
      REAL k_z1  = sqrt(k_rho / k_p);
      REAL k_z2  = k_z1 * k_ux;
      REAL k_z3  = k_z1 * k_uy;
      REAL k_z4  = k_z1 * k_uz;
      REAL k_z5  = k_z1 * k_p;

      REAL z1     = (REAL)(0.5) * (m_z1 + k_z1);
      REAL z1_log = logmean(m_z1, k_z1);
      REAL z2     = (REAL)(0.5) * (m_z2 + k_z2);
      REAL z3     = (REAL)(0.5) * (m_z3 + k_z3);
      REAL z4     = (REAL)(0.5) * (m_z4 + k_z4);
      REAL z5     = (REAL)(0.5) * (m_z5 + k_z5);
      REAL z5_log = logmean(m_z5, k_z5);

      REAL rho = z1 * z5_log;
      REAL ux = z2 / z1;
      REAL uy = z3 / z1;
      REAL uz = z4 / z1;
      REAL p1 = z5 / z1;
      REAL p2 = ((REAL)(GAMMA+1)*z5_log/z1_log + (REAL)(GAMMA-1)*z5/z1) / (REAL)(2*GAMMA);
      REAL h  = (REAL)(GAMMA/(GAMMA-1)) * p2/rho + (REAL)(0.5) * (ux*ux + uy*uy + uz*uz);

      REAL f_rho   = rho*ux;
      REAL f_rhoux = ux*f_rho + p1;
      REAL f_rhouy = uy*f_rho;
      REAL f_rhouz = uz*f_rho;
      REAL f_E     =  h*f_rho;

      ext_num_flux[Field_rho]    = - f_rho;
      ext_num_flux[Field_rho_ux] = - f_rhoux;
      ext_num_flux[Field_rho_uy] = - f_rhouy;
      ext_num_flux[Field_rho_uz] = - f_rhouz;
      ext_num_flux[Field_E]      = - f_E;
  }

  inline void compute_ext_num_flux_y(REAL const* uk, REAL const* um, REAL* ext_num_flux) {

      REAL m_rho = um[Field_rho];
      REAL m_E   = um[Field_E];
      REAL m_ux  = um[Field_ux];
      REAL m_uy  = um[Field_uy];
      REAL m_uz  = um[Field_uz];
      REAL m_p   = um[Field_p];
      REAL m_z1  = sqrt(m_rho / m_p);
      REAL m_z2  = m_z1 * m_ux;
      REAL m_z3  = m_z1 * m_uy;
      REAL m_z4  = m_z1 * m_uz;
      REAL m_z5  = m_z1 * m_p;

      REAL k_rho = uk[Field_rho];
      REAL k_E   = uk[Field_E];
      REAL k_ux  = uk[Field_ux];
      REAL k_uy  = uk[Field_uy];
      REAL k_uz  = uk[Field_uz];
      REAL k_p   = uk[Field_p];
      REAL k_z1  = sqrt(k_rho / k_p);
      REAL k_z2  = k_z1 * k_ux;
      REAL k_z3  = k_z1 * k_uy;
      REAL k_z4  = k_z1 * k_uz;
      REAL k_z5  = k_z1 * k_p;

      REAL z1     = (REAL)(0.5) * (m_z1 + k_z1);
      REAL z1_log = logmean(m_z1, k_z1);
      REAL z2     = (REAL)(0.5) * (m_z2 + k_z2);
      REAL z3     = (REAL)(0.5) * (m_z3 + k_z3);
      REAL z4     = (REAL)(0.5) * (m_z4 + k_z4);
      REAL z5     = (REAL)(0.5) * (m_z5 + k_z5);
      REAL z5_log = logmean(m_z5, k_z5);

      REAL rho = z1 * z5_log;
      REAL ux = z2 / z1;
      REAL uy = z3 / z1;
      REAL uz = z4 / z1;
      REAL p1 = z5 / z1;
      REAL p2 = ((REAL)(GAMMA+1)*z5_log/z1_log + (REAL)(GAMMA-1)*z5/z1) / (REAL)(2*GAMMA);
      REAL h  = (REAL)(GAMMA/(GAMMA-1)) * p2/rho + (REAL)(0.5) * (ux*ux + uy*uy + uz*uz);

      REAL f_rho   = rho*uy;
      REAL f_rhoux = ux*f_rho;
      REAL f_rhouy = uy*f_rho + p1;
      REAL f_rhouz = uz*f_rho;
      REAL f_E     =  h*f_rho;

      ext_num_flux[Field_rho]    = - f_rho;
      ext_num_flux[Field_rho_ux] = - f_rhoux;
      ext_num_flux[Field_rho_uy] = - f_rhouy;
      ext_num_flux[Field_rho_uz] = - f_rhouz;
      ext_num_flux[Field_E]      = - f_E;
  }

  inline void compute_ext_num_flux_z(REAL const* uk, REAL const* um, REAL* ext_num_flux) {

      REAL m_rho = um[Field_rho];
      REAL m_E   = um[Field_E];
      REAL m_ux  = um[Field_ux];
      REAL m_uy  = um[Field_uy];
      REAL m_uz  = um[Field_uz];
      REAL m_p   = um[Field_p];
      REAL m_z1  = sqrt(m_rho / m_p);
      REAL m_z2  = m_z1 * m_ux;
      REAL m_z3  = m_z1 * m_uy;
      REAL m_z4  = m_z1 * m_uz;
      REAL m_z5  = m_z1 * m_p;

      REAL k_rho = uk[Field_rho];
      REAL k_E   = uk[Field_E];
      REAL k_ux  = uk[Field_ux];
      REAL k_uy  = uk[Field_uy];
      REAL k_uz  = uk[Field_uz];
      REAL k_p   = uk[Field_p];
      REAL k_z1  = sqrt(k_rho / k_p);
      REAL k_z2  = k_z1 * k_ux;
      REAL k_z3  = k_z1 * k_uy;
      REAL k_z4  = k_z1 * k_uz;
      REAL k_z5  = k_z1 * k_p;

      REAL z1     = (REAL)(0.5) * (m_z1 + k_z1);
      REAL z1_log = logmean(m_z1, k_z1);
      REAL z2     = (REAL)(0.5) * (m_z2 + k_z2);
      REAL z3     = (REAL)(0.5) * (m_z3 + k_z3);
      REAL z4     = (REAL)(0.5) * (m_z4 + k_z4);
      REAL z5     = (REAL)(0.5) * (m_z5 + k_z5);
      REAL z5_log = logmean(m_z5, k_z5);

      REAL rho = z1 * z5_log;
      REAL ux = z2 / z1;
      REAL uy = z3 / z1;
      REAL uz = z4 / z1;
      REAL p1 = z5 / z1;
      REAL p2 = ((REAL)(GAMMA+1)*z5_log/z1_log + (REAL)(GAMMA-1)*z5/z1) / (REAL)(2*GAMMA);
      REAL h  = (REAL)(GAMMA/(GAMMA-1)) * p2/rho + (REAL)(0.5) * (ux*ux + uy*uy + uz*uz);

      REAL f_rho   = rho*uz;
      REAL f_rhoux = ux*f_rho;
      REAL f_rhouy = uy*f_rho;
      REAL f_rhouz = uz*f_rho + p1;
      REAL f_E     =  h*f_rho;

      ext_num_flux[Field_rho]    = - f_rho;
      ext_num_flux[Field_rho_ux] = - f_rhoux;
      ext_num_flux[Field_rho_uy] = - f_rhouy;
      ext_num_flux[Field_rho_uz] = - f_rhouz;
      ext_num_flux[Field_E]      = - f_E;
  }

#else

  #error "Error in ideal_gas_Euler.cl: No conservative discretisation specified!"

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


#if defined USE_BOUNDARY_FLUX_Suliciu
/* Suliciu relaxation solver (numerical flux) of
@book{bouchut2004nonlinear,
  title={Nonlinear Stability of Finite Volume Methods for Hyperbolic Conservation
         Laws and Well{\hyphen}Balanced Schemes for Sources},
  author={Bouchut, Fran{\c{c}}ois},
  year={2004},
  publisher={Birkh{\"a}user Verlag},
  address={Basel},
  doi={10.1007/b93802}
}
*/

  inline void compute_boundary_num_flux_x(REAL const* u_l, REAL const* u_r, REAL* num_flux) {

      // "left" state
      REAL l_rho = u_l[Field_rho];
      REAL l_E   = u_l[Field_E];
      REAL l_ux  = u_l[Field_ux];
      REAL l_uy  = u_l[Field_uy];
      REAL l_uz  = u_l[Field_uz];
      REAL l_p   = u_l[Field_p];
      REAL l_c = sqrt(GAMMA * l_p / l_rho);
      REAL l_eps = l_p / ((GAMMA-1) * l_rho);

      // "right" state
      REAL r_rho = u_r[Field_rho];
      REAL r_E   = u_r[Field_E];
      REAL r_ux  = u_r[Field_ux];
      REAL r_uy  = u_r[Field_uy];
      REAL r_uz  = u_r[Field_uz];
      REAL r_p   = u_r[Field_p];
      REAL r_c = sqrt(GAMMA * r_p / r_rho);
      REAL r_eps = r_p / ((GAMMA-1) * r_rho);

      REAL alpha = 0.5 * (GAMMA+1);

      // compute speeds
      REAL l_c_rho = 0;
      REAL r_c_rho = 0;
      if ((l_p <= r_p) && (0 < r_p)) {
        l_c_rho = l_c + alpha * fmax((REAL)(0), (r_p - l_p) / (r_rho * r_c) + l_ux - r_ux);
        r_c_rho = r_c + alpha * fmax((REAL)(0), (l_p - r_p) / (l_c_rho * l_rho) + l_ux - r_ux);
      }
      else if ((r_p <= l_p) && (0 < l_p)) {
        r_c_rho = r_c + alpha * fmax((REAL)(0), (l_p - r_p) / (l_rho * l_c) + l_ux - r_ux);
        l_c_rho = l_c + alpha * fmax((REAL)(0), (r_p - l_p) / (r_c_rho * r_rho) + l_ux - r_ux);
      }

      // compute intermediate values
      l_c = l_c_rho * l_rho;
      r_c = r_c_rho * r_rho;
      REAL s_ux = ( l_c*l_ux + r_c*r_ux + l_p - r_p) / (l_c + r_c);
      // = ifelse(isnan(... ?
      REAL s_p = (r_c*l_p + l_c*r_p - l_c*r_c*(r_ux - l_ux)) / (l_c + r_c);
      // = ifelse(isnan(... ?
      REAL ls_rho = (REAL)(1) / ( 1/l_rho + (r_c*(r_ux-l_ux) + l_p - r_p) / (l_c*(l_c+r_c)) );
      // = ifelse(isnan(... ?
      REAL rs_rho = (REAL)(1) / ( 1/r_rho + (l_c*(r_ux-l_ux) + r_p - l_p) / (r_c*(l_c+r_c)) );
      // = ifelse(isnan(... ?
      REAL ls_eps = l_eps + (s_p*s_p - l_p*l_p) / (2*l_c*l_c);
      REAL rs_eps = r_eps + (s_p*s_p - r_p*r_p) / (2*r_c*r_c);

      // compute fluxes
      // if (...) {...} else if (...) {...} etc.?
      REAL f_rho   = (0 <= l_ux-l_c_rho) *
                      (l_rho * l_ux)
                   + (l_ux-l_c_rho < 0) * (0 <= s_ux) *
                      (ls_rho * s_ux)
                   + (l_ux-l_c_rho < 0) *(s_ux < 0) * (0 <= r_ux+r_c_rho) *
                      (rs_rho * s_ux)
                   + (l_ux-l_c_rho < 0) *(s_ux < 0) * (r_ux+r_c_rho < 0) *
                      (r_rho * r_ux);
      REAL f_rhoux = (0 <= l_ux-l_c_rho) *
                      (l_rho * l_ux * l_ux + l_p)
                   + (l_ux-l_c_rho < 0) * (0 <= s_ux) *
                      (ls_rho * s_ux * s_ux + s_p)
                   + (l_ux-l_c_rho < 0) *(s_ux < 0) * (0 <= r_ux+r_c_rho) *
                      (rs_rho * s_ux * s_ux + s_p)
                   + (l_ux-l_c_rho < 0) *(s_ux < 0) * (r_ux+r_c_rho < 0) *
                      (r_rho * r_ux * r_ux + r_p);
      REAL f_rhouy = (0 <= l_ux-l_c_rho) *
                      (l_rho * l_ux * l_uy)
                   + (l_ux-l_c_rho < 0) * (0 <= s_ux) *
                      (ls_rho * s_ux * l_uy)
                   + (l_ux-l_c_rho < 0) *(s_ux < 0) * (0 <= r_ux+r_c_rho) *
                      (rs_rho * s_ux * r_uy)
                   + (l_ux-l_c_rho < 0) *(s_ux < 0) * (r_ux+r_c_rho < 0) *
                      (r_rho * r_ux * r_uy);
      REAL f_rhouz = (0 <= l_ux-l_c_rho) *
                      (l_rho * l_ux * l_uz)
                   + (l_ux-l_c_rho < 0) * (0 <= s_ux) *
                      (ls_rho * s_ux * l_uz)
                   + (l_ux-l_c_rho < 0) *(s_ux < 0) * (0 <= r_ux+r_c_rho) *
                      (rs_rho * s_ux * r_uz)
                   + (l_ux-l_c_rho < 0) *(s_ux < 0) * (r_ux+r_c_rho < 0) *
                      (r_rho * r_ux * r_uz);
      REAL f_E     = (0 <= l_ux-l_c_rho) *
                      (0.5 * l_rho * (l_ux*l_ux + l_uy*l_uy + l_uz*l_uz) + l_rho*l_eps + l_p) * l_ux
                   + (l_ux-l_c_rho < 0) * (0 <= s_ux) *
                      (0.5 * ls_rho * (s_ux*s_ux + l_uy*l_uy + l_uz*l_uz) + ls_rho*ls_eps + s_p) * s_ux
                   + (l_ux-l_c_rho < 0) *(s_ux < 0) * (0 <= r_ux+r_c_rho) *
                      (0.5 * ls_rho * (s_ux*s_ux + r_uy*r_uy + r_uz*r_uz) + ls_rho*ls_eps + s_p) * s_ux
                   + (l_ux-l_c_rho < 0) *(s_ux < 0) * (r_ux+r_c_rho < 0) *
                      (0.5 * r_rho * (r_ux*r_ux + r_uy*r_uy + r_uz*r_uz) + r_rho*r_eps + r_p) * r_ux;

      num_flux[Field_rho]    = f_rho;
      num_flux[Field_rho_ux] = f_rhoux;
      num_flux[Field_rho_uy] = f_rhouy;
      num_flux[Field_rho_uz] = f_rhouz;
      num_flux[Field_E]      = f_E;
  }

  inline void compute_boundary_num_flux_y(REAL const* u_l, REAL const* u_r, REAL* num_flux) {

      // "left" state
      REAL l_rho = u_l[Field_rho];
      REAL l_E   = u_l[Field_E];
      REAL l_ux  = u_l[Field_ux];
      REAL l_uy  = u_l[Field_uy];
      REAL l_uz  = u_l[Field_uz];
      REAL l_p   = u_l[Field_p];
      REAL l_c = sqrt(GAMMA * l_p / l_rho);
      REAL l_eps = l_p / ((GAMMA-1) * l_rho);

      // "right" state
      REAL r_rho = u_r[Field_rho];
      REAL r_E   = u_r[Field_E];
      REAL r_ux  = u_r[Field_ux];
      REAL r_uy  = u_r[Field_uy];
      REAL r_uz  = u_r[Field_uz];
      REAL r_p   = u_r[Field_p];
      REAL r_c = sqrt(GAMMA * r_p / r_rho);
      REAL r_eps = r_p / ((GAMMA-1) * r_rho);

      REAL alpha = 0.5 * (GAMMA+1);

      // compute speeds
      REAL l_c_rho = 0;
      REAL r_c_rho = 0;
      if ((l_p <= r_p) && (0 < r_p)) {
        l_c_rho = l_c + alpha * fmax((REAL)(0), (r_p - l_p) / (r_rho * r_c) + l_uy - r_uy);
        r_c_rho = r_c + alpha * fmax((REAL)(0), (l_p - r_p) / (l_c_rho * l_rho) + l_uy - r_uy);
      }
      else if ((r_p <= l_p) && (0 < l_p)) {
        r_c_rho = r_c + alpha * fmax((REAL)(0), (l_p - r_p) / (l_rho * l_c) + l_uy - r_uy);
        l_c_rho = l_c + alpha * fmax((REAL)(0), (r_p - l_p) / (r_c_rho * r_rho) + l_uy - r_uy);
      }

      // compute intermediate values
      l_c = l_c_rho * l_rho;
      r_c = r_c_rho * r_rho;
      REAL s_uy = ( l_c*l_uy + r_c*r_uy + l_p - r_p) / (l_c + r_c);
      // = ifelse(isnan(... ?
      REAL s_p = (r_c*l_p + l_c*r_p - l_c*r_c*(r_uy - l_uy)) / (l_c + r_c);
      //= ifelse(isnan(... ?
      REAL ls_rho = (REAL)(1) / ( 1/l_rho + (r_c*(r_uy-l_uy) + l_p - r_p) / (l_c*(l_c+r_c)) );
      // = ifelse(isnan(... ?
      REAL rs_rho = (REAL)(1) / ( 1/r_rho + (l_c*(r_uy-l_uy) + r_p - l_p) / (r_c*(l_c+r_c)) );
      // = ifelse(isnan(... ?
      REAL ls_eps = l_eps + (s_p*s_p - l_p*l_p) / (2*l_c*l_c);
      REAL rs_eps = r_eps + (s_p*s_p - r_p*r_p) / (2*r_c*r_c);

      // compute fluxes
      // if (...) {...} else if (...) {...} etc.?
      REAL f_rho   = (0 <= l_uy-l_c_rho) *
                      (l_rho * l_uy)
                   + (l_uy-l_c_rho < 0) * (0 <= s_uy) *
                      (ls_rho * s_uy)
                   + (l_uy-l_c_rho < 0) *(s_uy < 0) * (0 <= r_uy+r_c_rho) *
                      (rs_rho * s_uy)
                   + (l_uy-l_c_rho < 0) *(s_uy < 0) * (r_uy+r_c_rho < 0) *
                      (r_rho * r_uy);
      REAL f_rhoux = (0 <= l_uy-l_c_rho) *
                      (l_rho * l_uy * l_ux)
                   + (l_uy-l_c_rho < 0) * (0 <= s_uy) *
                      (ls_rho * s_uy * l_ux)
                   + (l_uy-l_c_rho < 0) *(s_uy < 0) * (0 <= r_uy+r_c_rho) *
                      (rs_rho * s_uy * r_ux)
                   + (l_uy-l_c_rho < 0) *(s_uy < 0) * (r_uy+r_c_rho < 0) *
                      (r_rho * r_uy * r_ux);
      REAL f_rhouy = (0 <= l_uy-l_c_rho) *
                      (l_rho * l_uy * l_uy + l_p)
                   + (l_uy-l_c_rho < 0) * (0 <= s_uy) *
                      (ls_rho * s_uy * l_uy + s_p)
                   + (l_uy-l_c_rho < 0) *(s_uy < 0) * (0 <= r_uy+r_c_rho) *
                      (rs_rho * s_uy * r_uy + s_p)
                   + (l_uy-l_c_rho < 0) *(s_uy < 0) * (r_uy+r_c_rho < 0) *
                      (r_rho * r_uy * r_uy + r_p);
      REAL f_rhouz = (0 <= l_uy-l_c_rho) *
                      (l_rho * l_uy * l_uz)
                   + (l_uy-l_c_rho < 0) * (0 <= s_uy) *
                      (ls_rho * s_uy * l_uz)
                   + (l_uy-l_c_rho < 0) *(s_uy < 0) * (0 <= r_uy+r_c_rho) *
                      (rs_rho * s_uy * r_uz)
                   + (l_uy-l_c_rho < 0) *(s_uy < 0) * (r_uy+r_c_rho < 0) *
                      (r_rho * r_uy * r_uz);
      REAL f_E     = (0 <= l_uy-l_c_rho) *
                      (0.5 * l_rho * (l_ux*l_ux + l_uy*l_uy + l_uz*l_uz) + l_rho*l_eps + l_p) * l_uy
                   + (l_uy-l_c_rho < 0) * (0 <= s_uy) *
                      (0.5 * ls_rho * (l_ux*l_ux + s_uy*s_uy + l_uz*l_uz) + ls_rho*ls_eps + s_p) * s_uy
                   + (l_uy-l_c_rho < 0) *(s_uy < 0) * (0 <= r_uy+r_c_rho) *
                      (0.5 * ls_rho * (r_ux*r_ux + s_uy*s_uy + r_uz*r_uz) + ls_rho*ls_eps + s_p) * s_uy
                   + (l_uy-l_c_rho < 0) *(s_uy < 0) * (r_uy+r_c_rho < 0) *
                      (0.5 * r_rho * (r_ux*r_ux + r_uy*r_uy + r_uz*r_uz) + r_rho*r_eps + r_p) * r_uy;

      num_flux[Field_rho]    = f_rho;
      num_flux[Field_rho_ux] = f_rhoux;
      num_flux[Field_rho_uy] = f_rhouy;
      num_flux[Field_rho_uz] = f_rhouz;
      num_flux[Field_E]      = f_E;
  }

  inline void compute_boundary_num_flux_z(REAL const* u_l, REAL const* u_r, REAL* num_flux) {

      // "left" state
      REAL l_rho = u_l[Field_rho];
      REAL l_E   = u_l[Field_E];
      REAL l_ux  = u_l[Field_ux];
      REAL l_uy  = u_l[Field_uy];
      REAL l_uz  = u_l[Field_uz];
      REAL l_p   = u_l[Field_p];
      REAL l_c = sqrt(GAMMA * l_p / l_rho);
      REAL l_eps = l_p / ((GAMMA-1) * l_rho);

      // "right" state
      REAL r_rho = u_r[Field_rho];
      REAL r_E   = u_r[Field_E];
      REAL r_ux  = u_r[Field_ux];
      REAL r_uy  = u_r[Field_uy];
      REAL r_uz  = u_r[Field_uz];
      REAL r_p   = u_r[Field_p];
      REAL r_c = sqrt(GAMMA * r_p / r_rho);
      REAL r_eps = r_p / ((GAMMA-1) * r_rho);

      REAL alpha = 0.5 * (GAMMA+1);

      // compute speeds
      REAL l_c_rho = 0;
      REAL r_c_rho = 0;
      if ((l_p <= r_p) && (0 < r_p)) {
        l_c_rho = l_c + alpha * fmax((REAL)(0), (r_p - l_p) / (r_rho * r_c) + l_uz - r_uz);
        r_c_rho = r_c + alpha * fmax((REAL)(0), (l_p - r_p) / (l_c_rho * l_rho) + l_uz - r_uz);
      }
      else if ((r_p <= l_p) && (0 < l_p)) {
        r_c_rho = r_c + alpha * fmax((REAL)(0), (l_p - r_p) / (l_rho * l_c) + l_uz - r_uz);
        l_c_rho = l_c + alpha * fmax((REAL)(0), (r_p - l_p) / (r_c_rho * r_rho) + l_uz - r_uz);
      }

      // compute intermediate values
      l_c = l_c_rho * l_rho;
      r_c = r_c_rho * r_rho;
      REAL s_uz = ( l_c*l_uz + r_c*r_uz + l_p - r_p) / (l_c + r_c);
      // = ifelse(isnan(... ?
      REAL s_p = (r_c*l_p + l_c*r_p - l_c*r_c*(r_uz - l_uz)) / (l_c + r_c);
      //= ifelse(isnan(... ?
      REAL ls_rho = (REAL)(1) / ( 1/l_rho + (r_c*(r_uz-l_uz) + l_p - r_p) / (l_c*(l_c+r_c)) );
      // = ifelse(isnan(... ?
      REAL rs_rho = (REAL)(1) / ( 1/r_rho + (l_c*(r_uz-l_uz) + r_p - l_p) / (r_c*(l_c+r_c)) );
      // = ifelse(isnan(... ?
      REAL ls_eps = l_eps + (s_p*s_p - l_p*l_p) / (2*l_c*l_c);
      REAL rs_eps = r_eps + (s_p*s_p - r_p*r_p) / (2*r_c*r_c);

      // compute fluxes
      // if (...) {...} else if (...) {...} etc.?
      REAL f_rho   = (0 <= l_uz-l_c_rho) *
                      (l_rho * l_uz)
                   + (l_uz-l_c_rho < 0) * (0 <= s_uz) *
                      (ls_rho * s_uz)
                   + (l_uz-l_c_rho < 0) *(s_uz < 0) * (0 <= r_uz+r_c_rho) *
                      (rs_rho * s_uz)
                   + (l_uz-l_c_rho < 0) *(s_uz < 0) * (r_uz+r_c_rho < 0) *
                      (r_rho * r_uz);
      REAL f_rhoux = (0 <= l_uz-l_c_rho) *
                      (l_rho * l_uz * l_ux)
                   + (l_uz-l_c_rho < 0) * (0 <= s_uz) *
                      (ls_rho * s_uz * l_ux)
                   + (l_uz-l_c_rho < 0) *(s_uz < 0) * (0 <= r_uz+r_c_rho) *
                      (rs_rho * s_uz * r_ux)
                   + (l_uz-l_c_rho < 0) *(s_uz < 0) * (r_uz+r_c_rho < 0) *
                      (r_rho * r_uz * r_ux);
      REAL f_rhouy = (0 <= l_uz-l_c_rho) *
                      (l_rho * l_uz * l_uy)
                   + (l_uz-l_c_rho < 0) * (0 <= s_uz) *
                      (ls_rho * s_uz * l_uy)
                   + (l_uz-l_c_rho < 0) *(s_uz < 0) * (0 <= r_uz+r_c_rho) *
                      (rs_rho * s_uz * r_uy)
                   + (l_uz-l_c_rho < 0) *(s_uz < 0) * (r_uz+r_c_rho < 0) *
                      (r_rho * r_uz * r_uy);
      REAL f_rhouz = (0 <= l_uz-l_c_rho) *
                      (l_rho * l_uz * l_uz + l_p)
                   + (l_uz-l_c_rho < 0) * (0 <= s_uz) *
                      (ls_rho * s_uz * l_uz + s_p)
                   + (l_uz-l_c_rho < 0) *(s_uz < 0) * (0 <= r_uz+r_c_rho) *
                      (rs_rho * s_uz * r_uz + s_p)
                   + (l_uz-l_c_rho < 0) *(s_uz < 0) * (r_uz+r_c_rho < 0) *
                      (r_rho * r_uz * r_uz + r_p);
      REAL f_E     = (0 <= l_uz-l_c_rho) *
                      (0.5 * l_rho * (l_ux*l_ux + l_uy*l_uy + l_uz*l_uz) + l_rho*l_eps + l_p) * l_uz
                   + (l_uz-l_c_rho < 0) * (0 <= s_uz) *
                      (0.5 * ls_rho * (l_ux*l_ux + l_uy*l_uy + s_uz*s_uz) + ls_rho*ls_eps + s_p) * s_uz
                   + (l_uz-l_c_rho < 0) *(s_uz < 0) * (0 <= r_uz+r_c_rho) *
                      (0.5 * ls_rho * (r_ux*r_ux + r_uy*r_uy + s_uz*s_uz) + ls_rho*ls_eps + s_p) * s_uz
                   + (l_uz-l_c_rho < 0) *(s_uz < 0) * (r_uz+r_c_rho < 0) *
                      (0.5 * r_rho * (r_ux*r_ux + r_uy*r_uy + r_uz*r_uz) + r_rho*r_eps + r_p) * r_uz;

      num_flux[Field_rho]    = f_rho;
      num_flux[Field_rho_ux] = f_rhoux;
      num_flux[Field_rho_uy] = f_rhouy;
      num_flux[Field_rho_uz] = f_rhouz;
      num_flux[Field_E]      = f_E;
  }

#else

  #error "Error in ideal_gas_Euler.cl: No boundary numerical flux specified!"

#endif


/* choice of surface term implementation:
  - USE_SURFACE_TERMS_MULTIPLICATION_WITH_BOOLS
  - USE_SURFACE_TERMS_IF
*/
#define USE_SURFACE_TERMS_IF


inline void add_surface_terms(REAL time, uint ix, uint iy, uint iz, global REAL *u, REAL *du_dt) {

  // For periodic boundary conditions and a single block, no surface term has to be used.
#ifndef USE_PERIODIC

  REAL u_inner[NUM_TOTAL_VARS] = {0};
  REAL u_outer[NUM_TOTAL_VARS] = {0};
  REAL flux[NUM_CONSERVED_VARS] = {0};
  REAL num_flux[NUM_CONSERVED_VARS] = {0};

  // inner and outer values of the conserved and auxiliary quantities
  get_field(ix, iy, iz, 0, 0, 0, u, u_inner);

  REAL rho_bound = rho_boundary(ix, iy, iz, time);
  REAL4 u_bound = u_boundary(ix, iy, iz, time);
  REAL p_bound = p_boundary(ix, iy, iz, time);
  REAL E_bound = compute_energy(p_bound, rho_bound, u_bound.x, u_bound.y, u_bound.z);
  u_outer[Field_rho] = rho_bound;
  u_outer[Field_rho_ux] = rho_bound*u_bound.x;
  u_outer[Field_rho_uy] = rho_bound*u_bound.y;
  u_outer[Field_rho_uz] = rho_bound*u_bound.z;
  u_outer[Field_E] = E_bound;
  u_outer[Field_ux] = u_bound.x;
  u_outer[Field_uy] = u_bound.y;
  u_outer[Field_uz] = u_bound.z;
  u_outer[Field_p] = p_bound;


#if defined USE_SURFACE_TERMS_MULTIPLICATION_WITH_BOOLS
  /*
  - Possibility 1

    |                       Device, Version                      | Runtime |
    |:----------------------------------------------------------:|--------:|
    | OpenCL 1.2 CUDA 9.1.84, GeForce GTX 1070 Ti                |   4.4 s |
    | OpenCL 2.1 LINUX, Intel(R) Core(TM) i7-8700K CPU @ 3.70GHz |  19.2 s |
    | OpenCL 2.1, Intel(R) Gen9 HD Graphics NEO                  |  22.6 s |
  */

  // flux x
  compute_flux_x(u_inner, flux);
  // left
  compute_boundary_num_flux_x(u_outer, u_inner, num_flux);
  for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
    du_dt[i] += check_bound_l(ix, 1) * (REAL)(M_INV[0] / DX) * (num_flux[i] - flux[i]);
  }
  // right
  compute_boundary_num_flux_x(u_inner, u_outer, num_flux);
  for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
    du_dt[i] -= check_bound_xr(ix, 1) * (REAL)(M_INV[0] / DX) * (num_flux[i] - flux[i]);
  }

  // flux y
  compute_flux_y(u_inner, flux);
  // bottom
  compute_boundary_num_flux_y(u_outer, u_inner, num_flux);
  for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
    du_dt[i] += check_bound_l(iy, 1) * (REAL)(M_INV[0] / DY) * (num_flux[i] - flux[i]);
  }
  // top
  compute_boundary_num_flux_y(u_inner, u_outer, num_flux);
  for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
    du_dt[i] -= check_bound_yr(iy, 1) * (REAL)(M_INV[0] / DY) * (num_flux[i] - flux[i]);
  }

  // flux z
  compute_flux_z(u_inner, flux);
  // left
  compute_boundary_num_flux_z(u_outer, u_inner, num_flux);
  for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
    du_dt[i] += check_bound_l(iz, 1) * (REAL)(M_INV[0] / DZ) * (num_flux[i] - flux[i]);
  }
  // right
  compute_boundary_num_flux_z(u_inner, u_outer, num_flux);
  for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
    du_dt[i] -= check_bound_zr(iz, 1) * (REAL)(M_INV[0] / DZ) * (num_flux[i] - flux[i]);
  }

#elif defined USE_SURFACE_TERMS_IF
  /*
  - Possibility 2

    |                       Device, Version                      | Runtime |
    |:----------------------------------------------------------:|--------:|
    | OpenCL 1.2 CUDA 9.1.84, GeForce GTX 1070 Ti                |   2.6 s |
    | OpenCL 2.1 LINUX, Intel(R) Core(TM) i7-8700K CPU @ 3.70GHz |  11.2 s |
    | OpenCL 2.1, Intel(R) Gen9 HD Graphics NEO                  |  12.4 s |
  */

  // flux x
  compute_flux_x(u_inner, flux);
  // left
  if (check_bound_l(ix, 1)) {
    compute_boundary_num_flux_x(u_outer, u_inner, num_flux);
    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
      du_dt[i] += (REAL)(M_INV[0] / DX) * (num_flux[i] - flux[i]);
    }
  }
  // right
  if (check_bound_xr(ix, 1)) {
    compute_boundary_num_flux_x(u_inner, u_outer, num_flux);
    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
      du_dt[i] -= (REAL)(M_INV[0] / DX) * (num_flux[i] - flux[i]);
    }
  }

  // flux y
  compute_flux_y(u_inner, flux);
  // bottom
  if (check_bound_l(iy, 1)) {
    compute_boundary_num_flux_y(u_outer, u_inner, num_flux);
    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
      du_dt[i] +=  (REAL)(M_INV[0] / DY) * (num_flux[i] - flux[i]);
    }
  }
  // top
  if (check_bound_yr(iy, 1)) {
    compute_boundary_num_flux_y(u_inner, u_outer, num_flux);
    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
      du_dt[i] -= (REAL)(M_INV[0] / DY) * (num_flux[i] - flux[i]);
    }
  }

  // flux z
  compute_flux_z(u_inner, flux);
  // left
  if (check_bound_l(iz, 1)) {
    compute_boundary_num_flux_z(u_outer, u_inner, num_flux);
    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
      du_dt[i] += (REAL)(M_INV[0] / DZ) * (num_flux[i] - flux[i]);
    }
  }
  // right
  if (check_bound_zr(iz, 1)) {
    compute_boundary_num_flux_z(u_inner, u_outer, num_flux);
    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
      du_dt[i] -= (REAL)(M_INV[0] / DZ) * (num_flux[i] - flux[i]);
    }
  }

#endif // USE_SURFACE_TERMS_...

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

  /*
  // No analytical solution
  for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
    set_field_component(ix, iy, iz, i, u, (REAL)(0));
  }*/

  REAL rho_ana = rho_analytical(ix, iy, iz, time);
  REAL4 u_ana = u_analytical(ix, iy, iz, time);
  REAL p_ana = p_analytical(ix, iy, iz, time);
  REAL E_ana = compute_energy(p_ana, rho_ana, u_ana.x, u_ana.y, u_ana.z);

  set_field_component(ix, iy, iz, Field_rho, u, rho_ana);
  set_field_component(ix, iy, iz, Field_rho_ux, u, rho_ana*u_ana.x);
  set_field_component(ix, iy, iz, Field_rho_uy, u, rho_ana*u_ana.y);
  set_field_component(ix, iy, iz, Field_rho_uz, u, rho_ana*u_ana.z);
  set_field_component(ix, iy, iz, Field_E, u, E_ana);
}

#endif // IDEAL_GAS_EULER
