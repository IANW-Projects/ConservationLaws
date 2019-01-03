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
