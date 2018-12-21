// This project is licensed under the terms of the Creative Commons CC BY-NC-ND 4.0 license.

#ifndef IDEAL_GAS_EULER_DEFINES_H
#define IDEAL_GAS_EULER_DEFINES_H

// ideal gas constant
REAL constant GAMMA = (REAL)(1.4);


#define NUM_CONSERVED_VARS 5
#define NUM_AUXILIARY_VARS 4
#define NUM_TOTAL_VARS (NUM_CONSERVED_VARS + NUM_AUXILIARY_VARS)

/* available numerical volume fluxes:
  - CENTRAL
  - DucrosEtAl

TODO: Finite difference fluxes
  - Morinishi
  - Kennedy and Gruber
  - Pirozzoli
TODO: Entropy conservative fluxes
  - Ismail and Roe
  - Chandrashekar
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
