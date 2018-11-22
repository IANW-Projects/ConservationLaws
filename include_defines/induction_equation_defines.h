// This project is licensed under the terms of the Creative Commons CC BY-NC-ND 4.0 license.

#ifndef INDUCTION_EQUATION_DEFINES_H
#define INDUCTION_EQUATION_DEFINES_H


#define NUM_CONSERVED_VARS 3
#define NUM_AUXILIARY_VARS 3
#define NUM_TOTAL_VARS (NUM_CONSERVED_VARS + NUM_AUXILIARY_VARS)

#define USE_UIBJ_PRODUCT
#define USE_SOURCE_CENTRAL
#define USE_UJBI_CENTRAL

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


#endif // INDUCTION_EQUATION_DEFINES_H