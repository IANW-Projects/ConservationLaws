// The high order artificial dissipation, ported here from InductionEq
// TODO: implement per field dissipation factors

#ifndef HO_DISSIPATION_FACTOR
	define HO_DISSIPATION_FACTOR 0.01
#endif

// Auxillary function for high order dissipation
// Used to compute the intermediate result D*a of equation (?)
void inline diff_HOD_x(uint ix, uint iy, uint iz, int bound_x, int start, global REAL const *d_field, REAL *HOD) {

  REAL um[NUM_TOTAL_VARS] = {0.0};
  int i, j;

  for (i = 0; i < STENCIL_WIDTH_HOD; i++) {
    get_field(ix, iy, iz, i + start, 0, 0, d_field, um);

    for (j = 0; i < NUM_CONSERVED_VARS; i++) {
    	HOD[j] = HOD[j]
        + M_INV_X[NUM_BOUNDS_HOD + bound_x]
        * D_HO[NUM_BOUNDS_HOD + bound_x][start + STENCIL_WIDTH_HOD - 1]
        * D_HO[NUM_BOUNDS_HOD][i] * pown(-1.0, (int)(ORDER / 2 + 1))
        * um[j];
    }
  }

}

void inline diff_HOD_y(uint ix, uint iy, uint iz, int bound_y, int start, global REAL const *d_field, REAL *HOD) {

  REAL um[NUM_TOTAL_VARS] = {0.0};
  int i, j;

  for (i = 0; i < STENCIL_WIDTH_HOD; i++) {
    get_field(ix, iy, iz, 0, i + start, 0, d_field, um);

    for (j = 0; j < NUM_CONSERVED_VARS; j++) {
        HOD[j] = HOD[j]
        + M_INV_Y[NUM_BOUNDS_HOD + bound_y]
        * D_HO[NUM_BOUNDS_HOD + bound_y][start + STENCIL_WIDTH_HOD - 1]
        * D_HO[NUM_BOUNDS_HOD][i] * pown(-1.0, (int)(ORDER / 2 + 1))
        * um[j];
    }   
  }
  
}

void inline diff_HOD_z(uint ix, uint iy, uint iz, int bound_z, int start, global REAL const *d_field, REAL *HOD) {

  REAL um[NUM_TOTAL_VARS] = {0.0};
  int i, j;

  for (i = 0; i < STENCIL_WIDTH_HOD; i++) {
    get_field(ix, iy, iz, 0,0, i + start, d_field, um);

    for (j = 0; j < NUM_CONSERVED_VARS; j++) {
        HOD[j] = HOD[j]
        + M_INV_Y[NUM_BOUNDS_HOD + bound_z]
        * D_HO[NUM_BOUNDS_HOD + bound_z][start + STENCIL_WIDTH_HOD - 1]
        * D_HO[NUM_BOUNDS_HOD][i] * pown(-1.0, (int)(ORDER / 2 + 1))
        * um[j];
    }
  }

}



// Computes the high order dissipation for fields in u at point (ix, iy, iz)
void inline high_order_dissipation(uint ix, uint iy, uint iz, global REAL const *u, REAL *HOD) {
  int i, j;
  uint idx = calc_idx(ix,iy,iz);

  int bound_x = get_bound_x(ix, NUM_BOUNDS_HOD);
  int bound_y = get_bound_y(iy, NUM_BOUNDS_HOD);
  int bound_z = get_bound_z(iz, NUM_BOUNDS_HOD);
  REAL HOD_X[NUM_TOTAL_VARS] = {0.0};
  REAL HOD_Y[NUM_TOTAL_VARS] = {0.0};
  REAL HOD_Z[NUM_TOTAL_VARS] = {0.0};

  for (i = -STENCIL_WIDTH_HOD + 1; i < 1; i++) {
    diff_HOD_x(ix, iy, iz, bound_x, i, u, HOD_X);
    diff_HOD_y(ix, iy, iz, bound_x, i, u, HOD_Y);
    diff_HOD_z(ix, iy, iz, bound_x, i, u, HOD_Z);
    for(j = 0; j < NUM_CONSERVED_VARS;j++) {
      HOD[j] = HOD[j] + (REAL)(HO_DISSIPATION_FACTOR) * HOD_X[j]
                      + (REAL)(HO_DISSIPATION_FACTOR) * HOD_Y[j]
                      + (REAL)(HO_DISSIPATION_FACTOR) * HOD_Z[j];
		}
        }

}

