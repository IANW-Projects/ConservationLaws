// This project is licensed under the terms of the Creative Commons CC BY-NC-ND 4.0 license.

#ifndef INDUCTION_EQUATION_HALL_DEFINES_H
#define INDUCTION_EQUATION_HALL_DEFINES_H


#define NUM_CONSERVED_VARS 3
#define NUM_AUXILIARY_VARS 7
#define NUM_TOTAL_VARS (NUM_CONSERVED_VARS + NUM_AUXILIARY_VARS)

#define USE_UIBJ_PRODUCT
#define USE_SOURCE_CENTRAL
#define USE_UJBI_CENTRAL

#define USE_HALL

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


REAL constant CONST_n1 = (REAL)(0.5773502691896258); // =1/sqrt(3)
REAL constant CONST_n2 = (REAL)(0.5773502691896258);
REAL constant CONST_n3 = (REAL)(0.5773502691896258);
REAL constant CONST_a = (REAL)(1.0);
REAL constant CONST_b = (REAL)(1.0);
REAL constant CONST_c = (REAL)(1.0);
REAL constant CONST_alpha = (REAL)(0.5);

#endif // INDUCTION_EQUATION_HALL_DEFINES_H
