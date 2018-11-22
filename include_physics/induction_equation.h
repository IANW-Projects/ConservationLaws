// This project is licensed under the terms of the Creative Commons CC BY-NC-ND 4.0 license.

#ifndef INDUCTION_EQUATION_H
#define INDUCTION_EQUATION_H

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
// Computation of auxiliary variables
//-------------------------------------------------------------------------------------------------------------------------------------------------------------

inline void compute_auxiliary_variables(REAL time, uint ix, uint iy, uint iz, global REAL *u) {

  REAL x = (REAL)XMIN + ix*(REAL)DX;
  REAL y = (REAL)YMIN + iy*(REAL)DY;

  set_field_component(ix, iy, iz, Field_ux, u, -y);
  set_field_component(ix, iy, iz, Field_uy, u, x);
  set_field_component(ix, iy, iz, Field_uz, u, 0.0);
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
// Extended numerical fluxes for the volume terms
//-------------------------------------------------------------------------------------------------------------------------------------------------------------

// 
/*
Extended numerical fluxes for '\partial_j(u_i B_j)' in space directions x, y, z.

Note: An aditional argument `uint dir` to access the vector via `Bk[dir]` instead
of `Bk.x` etc. can be used only for OpenCL version 2 and newer.
*/
void inline ext_num_flux_x_uiBj(REAL* uk, REAL* um, REAL* ext_num_flux) {

	#ifdef USE_UIBJ_CENTRAL
    ext_num_flux[Field_Bx] = ext_num_flux[Field_Bx] +  (REAL)(0.5) * (um[Field_ux]*um[Field_Bx] + uk[Field_ux]*uk[Field_Bx]);
    ext_num_flux[Field_By] = ext_num_flux[Field_By] +  (REAL)(0.5) * (um[Field_uy]*um[Field_Bx] + uk[Field_uy]*uk[Field_Bx]);
    ext_num_flux[Field_Bz] = ext_num_flux[Field_Bz] +  (REAL)(0.5) * (um[Field_uz]*um[Field_Bx] + uk[Field_uz]*uk[Field_Bx]);
	#elif defined USE_UIBJ_SPLIT
    ext_num_flux[Field_Bx] = ext_num_flux[Field_Bx] + (REAL)(0.25) * (um[Field_Bx] + uk[Field_Bx]) * (um[Field_ux] + uk[Field_ux]);
    ext_num_flux[Field_By] = ext_num_flux[Field_By] + (REAL)(0.25) * (um[Field_Bx] + uk[Field_Bx]) * (um[Field_uy] + uk[Field_uy]);
    ext_num_flux[Field_Bz] = ext_num_flux[Field_Bz] + (REAL)(0.25) * (um[Field_Bx] + uk[Field_Bx]) * (um[Field_uz] + uk[Field_uz]);
	#elif defined USE_UIBJ_PRODUCT
    ext_num_flux[Field_Bx] = ext_num_flux[Field_Bx] + (REAL)(0.5) * (um[Field_ux]*uk[Field_Bx] + uk[Field_ux]*um[Field_Bx]);
    ext_num_flux[Field_By] = ext_num_flux[Field_By] + (REAL)(0.5) * (um[Field_uy]*uk[Field_Bx] + uk[Field_uy]*um[Field_Bx]);
    ext_num_flux[Field_Bz] = ext_num_flux[Field_Bz] + (REAL)(0.5) * (um[Field_uz]*uk[Field_Bx] + uk[Field_uz]*um[Field_Bx]);
	#else
    #error "Error in induction_Eq_volume.cl: No discretization of 'D_j(u_i B_j)' specified!"
  #endif
}

void inline ext_num_flux_y_uiBj(REAL* uk, REAL* um, REAL* ext_num_flux) {

	#ifdef USE_UIBJ_CENTRAL
    ext_num_flux[Field_Bx] = ext_num_flux[Field_Bx] +  (REAL)(0.5) * (um[Field_ux]*um[Field_By] + uk[Field_ux]*uk[Field_By]);
    ext_num_flux[Field_By] = ext_num_flux[Field_By] +  (REAL)(0.5) * (um[Field_uy]*um[Field_By] + uk[Field_uy]*uk[Field_By]);
    ext_num_flux[Field_Bz] = ext_num_flux[Field_Bz] +  (REAL)(0.5) * (um[Field_uz]*um[Field_By] + uk[Field_uz]*uk[Field_By]);
	#elif defined USE_UIBJ_SPLIT
		ext_num_flux[Field_Bx] = ext_num_flux[Field_Bx] + (REAL)(0.25) * (um[Field_By] + uk[Field_By]) * (um[Field_ux] + uk[Field_ux]);
    ext_num_flux[Field_By] = ext_num_flux[Field_By] + (REAL)(0.25) * (um[Field_By] + uk[Field_By]) * (um[Field_uy] + uk[Field_uy]);
    ext_num_flux[Field_Bz] = ext_num_flux[Field_Bz] + (REAL)(0.25) * (um[Field_By] + uk[Field_By]) * (um[Field_uz] + uk[Field_uz]);
	#elif defined USE_UIBJ_PRODUCT
		ext_num_flux[Field_Bx] = ext_num_flux[Field_Bx] + (REAL)(0.5) * (um[Field_ux]*uk[Field_By] + uk[Field_ux]*um[Field_By]);
    ext_num_flux[Field_By] = ext_num_flux[Field_By] + (REAL)(0.5) * (um[Field_uy]*uk[Field_By] + uk[Field_uy]*um[Field_By]);
    ext_num_flux[Field_Bz] = ext_num_flux[Field_Bz] + (REAL)(0.5) * (um[Field_uz]*uk[Field_By] + uk[Field_uz]*um[Field_By]);
	#else
    #error "Error in induction_Eq_volume.cl: No discretization of 'D_j(u_i B_j)' specified!"
  #endif
}

void inline ext_num_flux_z_uiBj(REAL* uk, REAL* um, REAL* ext_num_flux) {

	#ifdef USE_UIBJ_CENTRAL
    ext_num_flux[Field_Bx] = ext_num_flux[Field_Bx] +  (REAL)(0.5) * (um[Field_ux]*um[Field_Bz] + uk[Field_ux]*uk[Field_Bz]);
    ext_num_flux[Field_By] = ext_num_flux[Field_By] +  (REAL)(0.5) * (um[Field_uy]*um[Field_Bz] + uk[Field_uy]*uk[Field_Bz]);
    ext_num_flux[Field_Bz] = ext_num_flux[Field_Bz] +  (REAL)(0.5) * (um[Field_uz]*um[Field_Bz] + uk[Field_uz]*uk[Field_Bz]);
	#elif defined USE_UIBJ_SPLIT
		ext_num_flux[Field_Bx] = ext_num_flux[Field_Bx] + (REAL)(0.25) * (um[Field_Bz] + uk[Field_Bz]) * (um[Field_ux] + uk[Field_ux]);
    ext_num_flux[Field_By] = ext_num_flux[Field_By] + (REAL)(0.25) * (um[Field_Bz] + uk[Field_Bz]) * (um[Field_uy] + uk[Field_uy]);
    ext_num_flux[Field_Bz] = ext_num_flux[Field_Bz] + (REAL)(0.25) * (um[Field_Bz] + uk[Field_Bz]) * (um[Field_uz] + uk[Field_uz]);
	#elif defined USE_UIBJ_PRODUCT
		ext_num_flux[Field_Bx] = ext_num_flux[Field_Bx] + (REAL)(0.5) * (um[Field_ux]*uk[Field_Bz] + uk[Field_ux]*um[Field_Bz]);
    ext_num_flux[Field_By] = ext_num_flux[Field_By] + (REAL)(0.5) * (um[Field_uy]*uk[Field_Bz] + uk[Field_uy]*um[Field_Bz]);
    ext_num_flux[Field_Bz] = ext_num_flux[Field_Bz] + (REAL)(0.5) * (um[Field_uz]*uk[Field_Bz] + uk[Field_uz]*um[Field_Bz]);
	#else
    #error "Error in induction_Eq_volume.cl: No discretization of 'D_j(u_i B_j)' specified!"
  #endif
}


/*
Extended numerical fluxes for 'u_i \partial_j B_j' in space directions x, y, z.

Note: An aditional argument `uint dir` to access the vector via `Bk[dir]` instead
of `Bk.x` etc. can be used only for OpenCL version 2 and newer.
*/
void inline ext_num_flux_x_source(REAL* uk, REAL* um, REAL* ext_num_flux) {

	#ifdef USE_SOURCE_ZERO
		ext_num_flux[Field_Bx] = ext_num_flux[Field_Bx] + 0;
    ext_num_flux[Field_By] = ext_num_flux[Field_By] + 0;
    ext_num_flux[Field_Bz] = ext_num_flux[Field_Bz] + 0;
	#elif defined USE_SOURCE_CENTRAL
		ext_num_flux[Field_Bx] = ext_num_flux[Field_Bx] + (REAL)(-0.5) * um[Field_ux] * (uk[Field_Bx] - um[Field_Bx]);
    ext_num_flux[Field_By] = ext_num_flux[Field_By] + (REAL)(-0.5) * um[Field_uy] * (uk[Field_Bx] - um[Field_Bx]);
    ext_num_flux[Field_Bz] = ext_num_flux[Field_Bz] + (REAL)(-0.5) * um[Field_uz] * (uk[Field_Bx] - um[Field_Bx]);
	#elif defined USE_SOURCE_SPLIT
		ext_num_flux[Field_Bx] = ext_num_flux[Field_Bx] + (REAL)(-0.25) * (um[Field_ux] + uk[Field_ux]) * (um[Field_Bx] + uk[Field_Bx]);
    ext_num_flux[Field_By] = ext_num_flux[Field_By] + (REAL)(-0.25) * (um[Field_uy] + uk[Field_uy]) * (um[Field_Bx] + uk[Field_Bx]);
    ext_num_flux[Field_Bz] = ext_num_flux[Field_Bz] + (REAL)(-0.25) * (um[Field_uz] + uk[Field_uz]) * (um[Field_Bx] + uk[Field_Bx]);
	#else
    #error "Error in induction_Eq_volume.cl: No discretization of 'u_i D_j B_j' specified!"
  #endif
}

REAL4 inline ext_num_flux_y_source(REAL* uk, REAL* um, REAL* ext_num_flux) {

	#ifdef USE_SOURCE_ZERO
		ext_num_flux[Field_Bx] = ext_num_flux[Field_Bx] + 0;
    ext_num_flux[Field_By] = ext_num_flux[Field_By] + 0;
    ext_num_flux[Field_Bz] = ext_num_flux[Field_Bz] + 0;
	#elif defined USE_SOURCE_CENTRAL
		ext_num_flux[Field_Bx] = ext_num_flux[Field_Bx] + (REAL)(-0.5) * um[Field_ux] * (uk[Field_By] - um[Field_By]);
    ext_num_flux[Field_By] = ext_num_flux[Field_By] + (REAL)(-0.5) * um[Field_uy] * (uk[Field_By] - um[Field_By]);
    ext_num_flux[Field_Bz] = ext_num_flux[Field_Bz] + (REAL)(-0.5) * um[Field_uz] * (uk[Field_By] - um[Field_By]);
	#elif defined USE_SOURCE_SPLIT
		ext_num_flux[Field_Bx] = ext_num_flux[Field_Bx] + (REAL)(-0.25) * (um[Field_ux] + uk[Field_ux]) * (um[Field_By] + uk[Field_By]);
    ext_num_flux[Field_By] = ext_num_flux[Field_By] + (REAL)(-0.25) * (um[Field_uy] + uk[Field_uy]) * (um[Field_By] + uk[Field_By]);
    ext_num_flux[Field_Bz] = ext_num_flux[Field_Bz] + (REAL)(-0.25) * (um[Field_uz] + uk[Field_uz]) * (um[Field_By] + uk[Field_By]);
	#else
    #error "Error in induction_Eq_volume.cl: No discretization of 'u_i D_j B_j' specified!"
  #endif
}

REAL4 inline ext_num_flux_z_source(REAL* uk, REAL* um, REAL* ext_num_flux) {

	#ifdef USE_SOURCE_ZERO
		ext_num_flux[Field_Bx] = ext_num_flux[Field_Bx] + 0;
    ext_num_flux[Field_By] = ext_num_flux[Field_By] + 0;
    ext_num_flux[Field_Bz] = ext_num_flux[Field_Bz] + 0;
	#elif defined USE_SOURCE_CENTRAL
		ext_num_flux[Field_Bx] = ext_num_flux[Field_Bx] + (REAL)(-0.5) * um[Field_ux] * (uk[Field_Bz] - um[Field_Bz]);
    ext_num_flux[Field_By] = ext_num_flux[Field_By] + (REAL)(-0.5) * um[Field_uy] * (uk[Field_Bz] - um[Field_Bz]);
    ext_num_flux[Field_Bz] = ext_num_flux[Field_Bz] + (REAL)(-0.5) * um[Field_uz] * (uk[Field_Bz] - um[Field_Bz]);
	#elif defined USE_SOURCE_SPLIT
		ext_num_flux[Field_Bx] = ext_num_flux[Field_Bx] + (REAL)(-0.25) * (um[Field_ux] + uk[Field_ux]) * (um[Field_Bz] + uk[Field_Bz]);
    ext_num_flux[Field_By] = ext_num_flux[Field_By] + (REAL)(-0.25) * (um[Field_uy] + uk[Field_uy]) * (um[Field_Bz] + uk[Field_Bz]);
    ext_num_flux[Field_Bz] = ext_num_flux[Field_Bz] + (REAL)(-0.25) * (um[Field_uz] + uk[Field_uz]) * (um[Field_Bz] + uk[Field_Bz]);
	#else
    #error "Error in induction_Eq_volume.cl: No discretization of 'u_i D_j B_j' specified!"
  #endif
}


/*
Extended numerical fluxes for '\partial_j(u_j B_i)' in space directions x, y, z.

Note: An aditional argument `uint dir` to access the vector via `Bk[dir]` instead
of `Bk.x` etc. can be used only for OpenCL version 2 and newer.
*/
void inline ext_num_flux_x_ujBi(REAL* uk, REAL* um, REAL* ext_num_flux) {

	#ifdef USE_UJBI_CENTRAL
  	ext_num_flux[Field_Bx] = ext_num_flux[Field_Bx] + (REAL)(-0.5) * (um[Field_ux]*um[Field_Bx] + uk[Field_ux]*uk[Field_Bx]);
    ext_num_flux[Field_By] = ext_num_flux[Field_By] + (REAL)(-0.5) * (um[Field_ux]*um[Field_By] + uk[Field_ux]*uk[Field_By]);
    ext_num_flux[Field_Bz] = ext_num_flux[Field_Bz] + (REAL)(-0.5) * (um[Field_ux]*um[Field_Bz] + uk[Field_ux]*uk[Field_Bz]);
	#elif defined USE_UJBI_SPLIT
    ext_num_flux[Field_Bx] = ext_num_flux[Field_Bx] + (REAL)(-0.25) * (um[Field_ux] + uk[Field_ux]) * (um[Field_Bx] + uk[Field_Bx]);
    ext_num_flux[Field_By] = ext_num_flux[Field_By] + (REAL)(-0.25) * (um[Field_ux] + uk[Field_ux]) * (um[Field_By] + uk[Field_By]);
    ext_num_flux[Field_Bz] = ext_num_flux[Field_Bz] + (REAL)(-0.25) * (um[Field_ux] + uk[Field_ux]) * (um[Field_Bz] + uk[Field_Bz]);
	#elif defined USE_UJBI_PRODUCT
    ext_num_flux[Field_Bx] = ext_num_flux[Field_Bx] + (REAL)(-0.5) * (um[Field_ux]*uk[Field_Bx] + uk[Field_ux]*um[Field_Bx]);
    ext_num_flux[Field_By] = ext_num_flux[Field_By] + (REAL)(-0.5) * (um[Field_ux]*uk[Field_By] + uk[Field_ux]*um[Field_By]);
    ext_num_flux[Field_Bz] = ext_num_flux[Field_Bz] + (REAL)(-0.5) * (um[Field_ux]*uk[Field_Bz] + uk[Field_ux]*um[Field_Bz]);
	#else
    #error "Error in induction_Eq_volume.cl: No discretization of 'D_j(u_j B_i)' specified!"
  #endif
}

void inline ext_num_flux_y_ujBi(REAL* uk, REAL* um, REAL* ext_num_flux) {

	#ifdef USE_UJBI_CENTRAL
  	ext_num_flux[Field_Bx] = ext_num_flux[Field_Bx] + (REAL)(-0.5) * (um[Field_uy]*um[Field_Bx] + uk[Field_uy]*uk[Field_Bx]);
    ext_num_flux[Field_By] = ext_num_flux[Field_By] + (REAL)(-0.5) * (um[Field_uy]*um[Field_By] + uk[Field_uy]*uk[Field_By]);
    ext_num_flux[Field_Bz] = ext_num_flux[Field_Bz] + (REAL)(-0.5) * (um[Field_uy]*um[Field_Bz] + uk[Field_uy]*uk[Field_Bz]);
	#elif defined USE_UJBI_SPLIT
    ext_num_flux[Field_Bx] = ext_num_flux[Field_Bx] + (REAL)(-0.25) * (um[Field_uy] + uk[Field_uy]) * (um[Field_Bx] + uk[Field_Bx]);
    ext_num_flux[Field_By] = ext_num_flux[Field_By] + (REAL)(-0.25) * (um[Field_uy] + uk[Field_uy]) * (um[Field_By] + uk[Field_By]);
    ext_num_flux[Field_Bz] = ext_num_flux[Field_Bz] + (REAL)(-0.25) * (um[Field_uy] + uk[Field_uy]) * (um[Field_Bz] + uk[Field_Bz]);
	#elif defined USE_UJBI_PRODUCT
    ext_num_flux[Field_Bx] = ext_num_flux[Field_Bx] + (REAL)(-0.5) * (um[Field_uy]*uk[Field_Bx] + uk[Field_uy]*um[Field_Bx]);
    ext_num_flux[Field_By] = ext_num_flux[Field_By] + (REAL)(-0.5) * (um[Field_uy]*uk[Field_By] + uk[Field_uy]*um[Field_By]);
    ext_num_flux[Field_Bz] = ext_num_flux[Field_Bz] + (REAL)(-0.5) * (um[Field_uy]*uk[Field_Bz] + uk[Field_uy]*um[Field_Bz]);
	#else
    #error "Error in induction_Eq_volume.cl: No discretization of 'D_j(u_j B_i)' specified!"
  #endif
}

void inline ext_num_flux_z_ujBi(REAL* uk, REAL* um, REAL* ext_num_flux) {

	#ifdef USE_UJBI_CENTRAL
  	ext_num_flux[Field_Bx] = ext_num_flux[Field_Bx] + (REAL)(-0.5) * (um[Field_uz]*um[Field_Bx] + uk[Field_uz]*uk[Field_Bx]);
    ext_num_flux[Field_By] = ext_num_flux[Field_By] + (REAL)(-0.5) * (um[Field_uz]*um[Field_By] + uk[Field_uz]*uk[Field_By]);
    ext_num_flux[Field_Bz] = ext_num_flux[Field_Bz] + (REAL)(-0.5) * (um[Field_uz]*um[Field_Bz] + uk[Field_uz]*uk[Field_Bz]);
	#elif defined USE_UJBI_SPLIT
    ext_num_flux[Field_Bx] = ext_num_flux[Field_Bx] + (REAL)(-0.25) * (um[Field_uz] + uk[Field_uz]) * (um[Field_Bx] + uk[Field_Bx]);
    ext_num_flux[Field_By] = ext_num_flux[Field_By] + (REAL)(-0.25) * (um[Field_uz] + uk[Field_uz]) * (um[Field_By] + uk[Field_By]);
    ext_num_flux[Field_Bz] = ext_num_flux[Field_Bz] + (REAL)(-0.25) * (um[Field_uz] + uk[Field_uz]) * (um[Field_Bz] + uk[Field_Bz]);
	#elif defined USE_UJBI_PRODUCT
    ext_num_flux[Field_Bx] = ext_num_flux[Field_Bx] + (REAL)(-0.5) * (um[Field_uz]*uk[Field_Bx] + uk[Field_uz]*um[Field_Bx]);
    ext_num_flux[Field_By] = ext_num_flux[Field_By] + (REAL)(-0.5) * (um[Field_uz]*uk[Field_By] + uk[Field_uz]*um[Field_By]);
    ext_num_flux[Field_Bz] = ext_num_flux[Field_Bz] + (REAL)(-0.5) * (um[Field_uz]*uk[Field_Bz] + uk[Field_uz]*um[Field_Bz]);
	#else
    #error "Error in induction_Eq_volume.cl: No discretization of 'D_j(u_j B_i)' specified!"
  #endif
}


/*
Extended numerical fluxes for the linear magnetic induction equation
'\partial_t B_i = \partial_j(u_i B_j - u_j B_i)'.

Note: An aditional argument `uint dir` to access the vector via `Bk[dir]` instead
of `Bk.x` etc. can be used only for OpenCL version 2 and newer.
*/
void inline compute_ext_num_flux_x(REAL* uk, REAL* um, REAL* ext_num_flux) {

	ext_num_flux_x_uiBj(uk, um, ext_num_flux);
  ext_num_flux_x_source(uk, um, ext_num_flux);
  ext_num_flux_x_ujBi(uk, um, ext_num_flux);
}

void inline compute_ext_num_flux_y(REAL* uk, REAL* um, REAL* ext_num_flux) {

	ext_num_flux_y_uiBj(uk, um, ext_num_flux);
  ext_num_flux_y_source(uk, um, ext_num_flux);
  ext_num_flux_y_ujBi(uk, um, ext_num_flux);
}

void inline compute_ext_num_flux_z(REAL* uk, REAL* um, REAL* ext_num_flux) {

	ext_num_flux_z_uiBj(uk, um, ext_num_flux);
  ext_num_flux_z_source(uk, um, ext_num_flux);
  ext_num_flux_z_ujBi(uk, um, ext_num_flux);
}


//--------------------------------------------------------------------------------------------------
// Functions for field initialization
//--------------------------------------------------------------------------------------------------

/*
Analytical solution of the magnetic field. It is used to calculate the numerical
error and related quantities.
*/
inline REAL4 rotation_2D_B(uint ix, uint iy, uint iz, REAL time) {

	REAL x = (REAL)XMIN + ix*(REAL)DX;
	REAL y = (REAL)YMIN + iy*(REAL)DY;
	REAL z = (REAL)ZMIN + iz*(REAL)DZ;

	REAL s = sin(time);
  REAL c = cos(time);

	REAL xx = c*x + s*y;
	REAL yy = -s*x + c*y;
	REAL zz = z;

	REAL fac = exp(-20*((xx-(REAL)(0.5))*(xx-(REAL)(0.5)) + yy*yy));
	REAL Bx0 = -4*yy*fac;
	REAL By0 = (4*xx-2)*fac;

	return (REAL4) {c*Bx0-s*By0, s*Bx0+c*By0, (REAL)(0), (REAL)(0)};
}


/*
Initial condition of the magnetic field.
*/
inline REAL4 b_init(uint ix, uint iy, uint iz) {

	return rotation_2D_B(ix, iy, iz, (REAL)0);
}

/*
Boundary condition of the magnetic field.
*/
inline REAL4 b_boundary(uint ix, uint iy, uint iz, REAL time) {

	return rotation_2D_B(ix, iy, iz, time);
}


inline void init_fields(uint ix, uint iy, uint iz, global REAL* u) {

  REAL4 B;
  B = b_init(ix, iy, iz);

  set_field_component(ix, iy, iz, Field_Bx, u, B.x);
  set_field_component(ix, iy, iz, Field_By, u, B.y);
  set_field_component(ix, iy, iz, Field_Bz, u, B.z);

}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
// Surface terms
//-------------------------------------------------------------------------------------------------------------------------------------------------------------


inline void add_surface_terms(REAL time, uint ix, uint iy, uint iz, global REAL *u, REAL *du_dt) {

  REAL um[NUM_TOTAL_VARS] = {0.0};;
  get_field(ix, iy, iz, 0, 0, 0, u, um);

// For periodic boundary conditions and a single block, no surface term has to be used.
#ifndef USE_PERIODIC
	REAL4 b_bound = b_boundary(ix, iy, iz, time);

  du_dt[Field_Bx] = du_dt[Field_Bx] 
                  + (REAL)(M_INV[0]/DX) * 
                  ((check_bound_l(ix,1) * (um[Field_ux] > 0)) * (REAL)(-1) + (check_bound_xr(ix,1) * (um[Field_ux] < 0)) * (REAL)(1)) * um[Field_ux] * (um[Field_Bx] - b_bound.x)
                  + (REAL)(M_INV[0]/DY) * 
                  ((check_bound_l(iy,1) * (um[Field_uy] > 0)) * (REAL)(-1) + (check_bound_yr(iy,1) * (um[Field_uy] < 0)) * (REAL)(1)) * um[Field_uy] * (um[Field_Bx] - b_bound.x)
                  + (REAL)(M_INV[0]/DZ) * 
                  ((check_bound_l(iz,1) * (um[Field_uz] > 0)) * (REAL)(-1) + (check_bound_zr(iz,1) * (um[Field_uz] < 0)) * (REAL)(1)) * um[Field_uz] * (um[Field_Bx] - b_bound.x);

  du_dt[Field_By] = du_dt[Field_By] 
                  + (REAL)(M_INV[0]/DX) * 
                  ((check_bound_l(ix,1) * (um[Field_ux] > 0)) * (REAL)(-1) + (check_bound_xr(ix,1) * (um[Field_ux] < 0)) * (REAL)(1)) * um[Field_ux] * (um[Field_By] - b_bound.y)
                  + (REAL)(M_INV[0]/DY) * 
                  ((check_bound_l(iy,1) * (um[Field_uy] > 0)) * (REAL)(-1) + (check_bound_yr(iy,1) * (um[Field_uy] < 0)) * (REAL)(1)) * um[Field_uy] * (um[Field_By] - b_bound.y)
                  + (REAL)(M_INV[0]/DZ) * 
                  ((check_bound_l(iz,1) * (um[Field_uz] > 0)) * (REAL)(-1) + (check_bound_zr(iz,1) * (um[Field_uz] < 0)) * (REAL)(1)) * um[Field_uz] * (um[Field_By] - b_bound.y);

  du_dt[Field_Bz] = du_dt[Field_Bz] 
                  + (REAL)(M_INV[0]/DX) * 
                  ((check_bound_l(ix,1) * (um[Field_ux] > 0)) * (REAL)(-1) + (check_bound_xr(ix,1) * (um[Field_ux] < 0)) * (REAL)(1)) * um[Field_ux] * (um[Field_Bz] - b_bound.z)
                  + (REAL)(M_INV[0]/DY) * 
                  ((check_bound_l(iy,1) * (um[Field_uy] > 0)) * (REAL)(-1) + (check_bound_yr(iy,1) * (um[Field_uy] < 0)) * (REAL)(1)) * um[Field_uy] * (um[Field_Bz] - b_bound.z)
                  + (REAL)(M_INV[0]/DZ) * 
                  ((check_bound_l(iz,1) * (um[Field_uz] > 0)) * (REAL)(-1) + (check_bound_zr(iz,1) * (um[Field_uz] < 0)) * (REAL)(1)) * um[Field_uz] * (um[Field_Bz] - b_bound.z);
#endif // USE_PERIODIC

}


#endif // INDUCTION_EQUATION_H
