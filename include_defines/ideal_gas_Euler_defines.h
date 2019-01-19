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
    - Ranocha (2018)
  - other fluxes
    - CENTRAL
    - DucrosEtAl (2000)
    - KennedyGruber (2008)
    - Morinishi (2010)
    - Pirozzoli (2011)
    - KuyaTotaniKawai_QCC (2018)

TODO: Finite difference fluxes
  - KravchenkoMoin (1997)
TODO: Entropy conservative fluxes
  - Variants of Hicken, Crean (2018)
*/
#define USE_FLUX_KennedyGruber


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
