// This project is licensed under the terms of the Creative Commons CC BY-NC-ND 4.0 license.

#ifndef IDEAL_MHD_DEFINES_H
#define IDEAL_MHD_DEFINES_H


#define NUM_CONSERVED_VARS 8
#define NUM_AUXILIARY_VARS 4
#define NUM_TOTAL_VARS (NUM_CONSERVED_VARS + NUM_AUXILIARY_VARS)

/* available numerical volume fluxes:
  - CENTRAL
  - SJOGREEN_YEE_DS2

TODO: SJOGREEN_YEE_DS1, SJOGREEN_YEE_DS3, SJOGREEN_YEE_DS4, cf.
@inproceedings{sjogreen2017skew,
  title={Skew-symmetric splitting and stability of high order central schemes},
  author={Sj{\"o}green, Bj{\"o}rn and Yee, Helen C and Kotov, Dmitry},
  booktitle={Journal of Physics: Conference Series},
  volume={837},
  number={1},
  pages={012019},
  year={2017},
  organization={IOP Publishing},
  doi={10.1088/1742-6596/837/1/012019}
}
TODO: Entropy conservative fluxes described in
@article{sjogreen2018high,
  title={High order entropy conservative central schemes for wide ranges of
         compressible gas dynamics and {MHD} flows},
  author={Sj{\"o}green, Bj{\"o}rn and Yee, HC},
  journal={Journal of Computational Physics},
  volume={364},
  pages={153--185},
  year={2018},
  publisher={Elsevier},
  doi={10.1016/j.jcp.2018.02.003}
}
*/
#define USE_FLUX_SJOGREEN_YEE_DS2
/* available source terms:
  - NONE
  - GODUNOV_CENTRAL

TODO: JANHUNEN_CENTRAL, BRACKBILL_BARNES_CENTRAL, cf.
@inproceedings{sjogreen2017skew,
  title={Skew-symmetric splitting and stability of high order central schemes},
  author={Sj{\"o}green, Bj{\"o}rn and Yee, Helen C and Kotov, Dmitry},
  booktitle={Journal of Physics: Conference Series},
  volume={837},
  number={1},
  pages={012019},
  year={2017},
  organization={IOP Publishing},
  doi={10.1088/1742-6596/837/1/012019}
}
*/
#define USE_SOURCE_GODUNOV_CENTRAL

// Boundary Fluxes: USE_BOUNDARY_FLUX_HLL
#define USE_BOUNDARY_FLUX_HLL

// Speeds: USE_HLL_SPEED_KUSANO USE_HLL_SPEED_BKW
#define USE_HLL_SPEED_BKW


enum Fields {
             // Density
             Field_rho,
             // Momentum
             Field_rho_ux,
             Field_rho_uy,
             Field_rho_uz,
             // Total Energy
             Field_E,
             // Magnetic Field
             Field_Bx,
             Field_By,
             Field_Bz,
             // Velocity
             Field_ux,
             Field_uy,
             Field_uz,
             // Pressure
             Field_p
            };


#endif // IDEAL_MHD_DEFINES_H
