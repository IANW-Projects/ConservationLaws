// This project is licensed under the terms of the Creative Commons CC BY-NC-ND 4.0 license.

#ifndef LINEAR_VARIABLE_ADVECTION_DEFINES_H
#define LINEAR_VARIABLE_ADVECTION_DEFINES_H


#define NUM_CONSERVED_VARS 1
#define NUM_AUXILIARY_VARS 3
#define NUM_TOTAL_VARS (NUM_CONSERVED_VARS + NUM_AUXILIARY_VARS)

// available numerical volume fluxes: central, split, product
#define USE_SPLIT_FLUX

enum Fields {
             // Unknown
             Field_u,
             // Velocity Field
             Field_ax,
             Field_ay,
             Field_az
            };


#endif // LINEAR_VARIABLE_ADVECTION_DEFINES_H
