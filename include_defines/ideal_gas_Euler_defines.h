// This project is licensed under the terms of the Creative Commons CC BY-NC-ND 4.0 license.

#ifndef IDEAL_GAS_EULER_DEFINES_H
#define IDEAL_GAS_EULER_DEFINES_H

// ideal gas constant
REAL constant GAMMA = (REAL)(1.4);


#define NUM_CONSERVED_VARS 5
#define NUM_AUXILIARY_VARS 4
#define NUM_TOTAL_VARS (NUM_CONSERVED_VARS + NUM_AUXILIARY_VARS)

/* available numerical volume fluxes:
  - entropy conservative fluxes
    - IsmailRoe (2009)
    - Chandrashekar (2013)
  - other fluxes
    - CENTRAL
    - DucrosEtAl (2000)
    - Morinishi (2010)
    - KuyaTotaniKawai_QCC (2018)

TODO: Finite difference fluxes
  - Kennedy and Gruber (2008)
  - Pirozzoli (2011)
TODO: Entropy conservative fluxes
  - Ranocha
*/
#define USE_FLUX_DucrosEtAl


enum Fields {
             // Density
             Field_rho,
             // Momentum
             Field_rho_ux,
             Field_rho_uy,
             Field_rho_uz,
             // Total Energy
             Field_E,
             // Velocity
             Field_ux,
             Field_uy,
             Field_uz,
             // Pressure
             Field_p
            };


#endif // IDEAL_GAS_EULER_DEFINES_H
