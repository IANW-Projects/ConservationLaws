// This project is licensed under the terms of the Creative Commons CC BY-NC-ND 4.0 license.

#ifndef INDUCTION_EQUATION_DEFINES_H
#define INDUCTION_EQUATION_DEFINES_H

#define USE_UIBJ_PRODUCT
#define USE_SOURCE_CENTRAL
#define USE_UJBI_CENTRAL

#define USE_HALL

#ifdef USE_HALL
    #define NUM_CONSERVED_VARS 3
    #define NUM_AUXILIARY_VARS 7
    #define NUM_TOTAL_VARS (NUM_CONSERVED_VARS + NUM_AUXILIARY_VARS)

    enum Fields {
             // Magnetic Field
             Field_Bx,
             Field_By,
             Field_Bz,
             // Velocity Field
             Field_ux,
             Field_uy,
             Field_uz,
             // Hall
             Field_curlB_rho_x,
             Field_curlB_rho_y,
             Field_curlB_rho_z,
             // Density
             Field_rho
            };
#else
    #define NUM_CONSERVED_VARS 3
    #define NUM_AUXILIARY_VARS 3
    #define NUM_TOTAL_VARS (NUM_CONSERVED_VARS + NUM_AUXILIARY_VARS)

    enum Fields {
                // Magnetic Field
                Field_Bx,
                Field_By,
                Field_Bz,
                // Velocity Field
                Field_ux,
                Field_uy,
                Field_uz
                };
#endif


#endif // INDUCTION_EQUATION_DEFINES_H