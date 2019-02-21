// This project is licensed under the terms of the Creative Commons CC BY-NC-ND 4.0 license.

// Contains kernel and auxillary kernel for time integration.


// Add the volume terms of an SBP SAT discretisation using extended numerical fluxes to du_dt.
inline void add_volume_terms(REAL time, uint ix, uint iy, uint iz, global REAL const* u, REAL *du_dt) {

  int bound_x = get_bound_x(ix, NUM_BOUNDS);
  int bound_y = get_bound_y(iy, NUM_BOUNDS);
  int bound_z = get_bound_z(iz, NUM_BOUNDS);

  // the variables at the current position
  REAL um[NUM_TOTAL_VARS] = {0.0};
  get_field(ix, iy, iz, 0, 0, 0, u, um);

  // the variables at remote positions used to compute the time derivative at (ix, iy, iz)
  REAL uk[NUM_TOTAL_VARS] = {0.0};

  // the time derivative of the conserved variables, initalised as zero
  REAL du[NUM_CONSERVED_VARS] = {0.0};

  // the extended numerical flux used to calculate the time derivative
  REAL ext_num_flux[NUM_CONSERVED_VARS] = {0.0};


  // x-direction
  for (uint j = 0; j < STENCIL_WIDTH; j++) {

    get_field(ix, iy, iz, (j - (STENCIL_WIDTH - 1)/2), 0, 0, u, uk);

    compute_ext_num_flux_x(uk, um, ext_num_flux);

    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
      du[i] = du[i] + SBP_diff[NUM_BOUNDS + bound_x][j] * ext_num_flux[i];
      ext_num_flux[i] = (REAL)(0);
    }
	}
  for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
    du_dt[i] = du_dt[i] + (REAL)(2.0/DX) * du[i];
    du[i] = (REAL)(0);
  }

  // y-direction
	for (uint j = 0; j < STENCIL_WIDTH; j++) {

    get_field(ix, iy, iz, 0, (j - (STENCIL_WIDTH - 1)/2), 0, u, uk);

    compute_ext_num_flux_y(uk, um, ext_num_flux);

    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
      du[i] = du[i] + SBP_diff[NUM_BOUNDS + bound_y][j] * ext_num_flux[i];
      ext_num_flux[i] = (REAL)(0);
    }
	}
  for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
    du_dt[i] = du_dt[i] + (REAL)(2.0/DY) * du[i];
    du[i] = (REAL)(0);
  }

  // z-direction
	for (uint j = 0; j < STENCIL_WIDTH; j++) {

  	get_field(ix, iy, iz, 0, 0, (j - (STENCIL_WIDTH - 1)/2), u, uk);

    compute_ext_num_flux_z(uk, um, ext_num_flux);

    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
      du[i] = du[i] + SBP_diff[NUM_BOUNDS + bound_z][j] * ext_num_flux[i];
      ext_num_flux[i] = (REAL)(0);
    }
	}
  for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
	  du_dt[i] = du_dt[i] + (REAL)(2.0/DZ) * du[i];
  }
}


// Compute the time derivative du_dt using an SBP SAT discretisation with extended numerical fluxes.
inline void compute_du_dt(REAL time, uint ix, uint iy, uint iz, global REAL const* u, REAL *du_dt) {

  for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
	  du_dt[i] = 0;
  }

  add_volume_terms(time, ix, iy, iz, u, du_dt);
  add_surface_terms(time, ix, iy, iz, u, du_dt);
}

//--------------------------------------------------------------------------------------------------
// Time update kernels
//--------------------------------------------------------------------------------------------------

/*
Update the global time variable as time -> time + DT. This kernel takes two arguments d_field_bi
and can be used with CarpenterKennedy2N54.
*/
kernel void calc_time_2_args(global REAL *u1, global REAL *u2, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    if (ix == 0 && iy == 0 && iz == 0) {
	      // time[0] is the time at the beginning of a time step
        // time[1] is the local time at one stage of the Runge-Kutta method
	      time[0] = time[0] + (REAL)(DT);
        time[1] = time[0];
	}
}


kernel void calc_time_3_args(global REAL *u1, global REAL *u2, global REAL *u3, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    if (ix == 0 && iy == 0 && iz == 0) {
	      // time[0] is the time at the beginning of a time step
        // time[1] is the local time at one stage of the Runge-Kutta method
	      time[0] = time[0] + (REAL)(DT);
        time[1] = time[0];
	}
}

//--------------------------------------------------------------------------------------------------
// Auxiliary time integrator kernels
//--------------------------------------------------------------------------------------------------

/*
Update the auxiliary variables. These kernels take two arguments ui and can be used with
CarpenterKennedy2N54.
*/
kernel void calc_auxiliary_vars_1_2_args(global REAL *u1, global REAL *u2, global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  compute_auxiliary_variables(time[1], ix, iy, iz, u1);
}

kernel void calc_auxiliary_vars_2_2_args(global REAL *u1, global REAL *u2, global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  compute_auxiliary_variables(time[1], ix, iy, iz, u2);
}

kernel void calc_auxiliary_vars_1_3_args(global REAL *u1, global REAL *u2, global REAL *u3, global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  compute_auxiliary_variables(time[1], ix, iy, iz, u1);
}

kernel void calc_auxiliary_vars_2_3_args(global REAL *u1, global REAL *u2, global REAL *u3, global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  compute_auxiliary_variables(time[1], ix, iy, iz, u2);
}

kernel void calc_auxiliary_vars_3_3_args(global REAL *u1, global REAL *u2, global REAL *u3, global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  compute_auxiliary_variables(time[1], ix, iy, iz, u3);
}


//--------------------------------------------------------------------------------------------------
// Time integrator kernels using 2 fields
//--------------------------------------------------------------------------------------------------


/*
	CarpenterKennedy2N54_X

Perform one time step using the five stage, fourth order, explicit Runge-Kutta method
CarpenterKennedy2N54. Starting with the old value `u1`, the new value after one time step is
obtained using the temporary array `u2`.
*/
kernel void CarpenterKennedy2N54_1a(global REAL *u1, global REAL *u2, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    REAL du[NUM_CONSERVED_VARS];
    compute_du_dt(time[1], ix, iy, iz, u1, du);

    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
        set_field_component(ix, iy, iz, i, u2, (REAL)(DT) * du[i]);
    }
}
kernel void CarpenterKennedy2N54_1b(global REAL *u1, global REAL *u2, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    REAL coef_B = (REAL)(1432997174477.0 / 9575080441755.0);
    REAL coef_C = (REAL)(1432997174477.0 / 9575080441755.0);

    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
        axpy_field_component(ix, iy, iz, i, (REAL)(1), u1, coef_B * get_field_component(ix, iy, iz, i, u2));
    }

    if (ix == 0 && iy == 0 && iz == 0) {
        // time[0] is the time at the beginning of a time step
        // time[1] is the local time at one stage of a Runge-Kutta method
        time[1] = time[0] + coef_C * (REAL)(DT);
    }
}

kernel void CarpenterKennedy2N54_2a(global REAL *u1, global REAL *u2, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    REAL du[NUM_CONSERVED_VARS];

    REAL coef_A = (REAL)(-567301805773.0 / 1357537059087.0);

    compute_du_dt(time[1], ix, iy, iz, u1, du);

    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
        axpy_field_component(ix, iy, iz, i, coef_A, u2, (REAL)(DT) * du[i]);
    }
}
kernel void CarpenterKennedy2N54_2b(global REAL *u1, global REAL *u2, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    REAL coef_B = (REAL)(5161836677717.0 / 13612068292357.0);
    REAL coef_C = (REAL)(2526269341429.0 / 6820363962896.0);

    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
        axpy_field_component(ix, iy, iz, i, (REAL)(1), u1, coef_B * get_field_component(ix, iy, iz, i, u2));
    }

    if (ix == 0 && iy == 0 && iz == 0) {
        // time[0] is the time at the beginning of a time step
        // time[1] is the local time at one stage of a Runge-Kutta method
        time[1] = time[0] + coef_C * (REAL)(DT);
    }
}

kernel void CarpenterKennedy2N54_3a(global REAL *u1, global REAL *u2, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    REAL du[NUM_CONSERVED_VARS] = {0.0};

    REAL coef_A = (REAL)(-2404267990393.0 / 2016746695238.0);

    compute_du_dt(time[1], ix, iy, iz, u1, du);

    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
        axpy_field_component(ix, iy, iz, i, coef_A, u2, (REAL)(DT) * du[i]);
    }
}
kernel void CarpenterKennedy2N54_3b(global REAL *u1, global REAL *u2, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    REAL coef_B = (REAL)(1720146321549.0 / 2090206949498.0);
    REAL coef_C = (REAL)(2006345519317.0 / 3224310063776.0);

    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
        axpy_field_component(ix, iy, iz, i, (REAL)(1), u1, coef_B * get_field_component(ix, iy, iz, i, u2));
    }

    if (ix == 0 && iy == 0 && iz == 0) {
        // time[0] is the time at the beginning of a time step
        // time[1] is the local time at one stage of a Runge-Kutta method
        time[1] = time[0] + coef_C * (REAL)(DT);
    }
}

kernel void CarpenterKennedy2N54_4a(global REAL *u1, global REAL *u2, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    REAL du[NUM_CONSERVED_VARS] = {0.0};

    REAL coef_A = (REAL)(-3550918686646.0 / 2091501179385.0);

    compute_du_dt(time[1], ix, iy, iz, u1, du);

    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
        axpy_field_component(ix, iy, iz, i, coef_A, u2, (REAL)(DT) * du[i]);
    }
}
kernel void CarpenterKennedy2N54_4b(global REAL *u1, global REAL *u2, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    REAL coef_B = (REAL)(3134564353537.0 / 4481467310338.0);
    REAL coef_C = (REAL)(2802321613138.0 / 2924317926251.0);

    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
        axpy_field_component(ix, iy, iz, i, (REAL)(1), u1, coef_B * get_field_component(ix, iy, iz, i, u2));
    }

    if (ix == 0 && iy == 0 && iz == 0) {
        // time[0] is the time at the beginning of a time step
        // time[1] is the local time at one stage of a Runge-Kutta method
        time[1] = time[0] + coef_C * (REAL)(DT);
    }
}

kernel void CarpenterKennedy2N54_5a(global REAL *u1, global REAL *u2, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    REAL du[NUM_CONSERVED_VARS] = {0.0};

    REAL coef_A = (REAL)(-1275806237668.0 / 842570457699.0);

    compute_du_dt(time[1], ix, iy, iz, u1, du);

    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
        axpy_field_component(ix, iy, iz, i, coef_A, u2, (REAL)(DT) * du[i]);
    }
    }
    kernel void CarpenterKennedy2N54_5b(global REAL *u1, global REAL *u2, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    REAL coef_B = (REAL)(2277821191437.0 / 14882151754819.0);

    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
        axpy_field_component(ix, iy, iz, i, (REAL)(1), u1, coef_B * get_field_component(ix, iy, iz, i, u2));
    }

    if (ix == 0 && iy == 0 && iz == 0) {
        // time[0] is the time at the beginning of a time step
        // time[1] is the local time at one stage of a Runge-Kutta method
        time[0] = time[0] + (REAL)(DT);
        time[1] = time[0];
    }
}

/*
    ToulorgeDesmet2N84F

Perform one time step using the eight stage, fourth order, explicit Runge-Kutta method
ToulorgeDesmet2N84F.
*/
kernel void ToulorgeDesmet2N84F_1a(global REAL *u1, global REAL *u2, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    REAL du[NUM_CONSERVED_VARS];
    compute_du_dt(time[1], ix, iy, iz, u1, du);

    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
        set_field_component(ix, iy, iz, i, u2, (REAL)(DT) * du[i]);
    }
}
kernel void ToulorgeDesmet2N84F_1b(global REAL *u1, global REAL *u2, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    REAL coef_B = (REAL)(0.08037936882736950);
    REAL coef_C = (REAL)(0.08037936882736950);

    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
        axpy_field_component(ix, iy, iz, i, (REAL)(1), u1, coef_B * get_field_component(ix, iy, iz, i, u2));
    }

    if (ix == 0 && iy == 0 && iz == 0) {
        // time[0] is the time at the beginning of a time step
        // time[1] is the local time at one stage of a Runge-Kutta method
        time[1] = time[0] + coef_C * (REAL)(DT);
	}
}

kernel void ToulorgeDesmet2N84F_2a(global REAL *u1, global REAL *u2, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    REAL coef_A = (REAL)(-0.5534431294501569);

    REAL du[NUM_CONSERVED_VARS];
    compute_du_dt(time[1], ix, iy, iz, u1, du);

    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
        axpy_field_component(ix, iy, iz, i, coef_A, u2, (REAL)(DT) * du[i]);
    }
}
kernel void ToulorgeDesmet2N84F_2b(global REAL *u1, global REAL *u2, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    REAL coef_B = (REAL)(0.5388497458569843);
    REAL coef_C = (REAL)(0.3210064250338430);

    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
        axpy_field_component(ix, iy, iz, i, (REAL)(1), u1, coef_B * get_field_component(ix, iy, iz, i, u2));
    }

    if (ix == 0 && iy == 0 && iz == 0) {
        // time[0] is the time at the beginning of a time step
        // time[1] is the local time at one stage of a Runge-Kutta method
        time[1] = time[0] + coef_C * (REAL)(DT);
	}
}

kernel void ToulorgeDesmet2N84F_3a(global REAL *u1, global REAL *u2, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    REAL coef_A = (REAL)(0.01065987570203490);

    REAL du[NUM_CONSERVED_VARS];
    compute_du_dt(time[1], ix, iy, iz, u1, du);

    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
        axpy_field_component(ix, iy, iz, i, coef_A, u2, (REAL)(DT) * du[i]);
    }
}
kernel void ToulorgeDesmet2N84F_3b(global REAL *u1, global REAL *u2, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    REAL coef_B = (REAL)(0.01974974409031960);
    REAL coef_C = (REAL)(0.3408501826604660);

    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
        axpy_field_component(ix, iy, iz, i, (REAL)(1), u1, coef_B * get_field_component(ix, iy, iz, i, u2));
    }

    if (ix == 0 && iy == 0 && iz == 0) {
        // time[0] is the time at the beginning of a time step
        // time[1] is the local time at one stage of a Runge-Kutta method
        time[1] = time[0] + coef_C * (REAL)(DT);
	}
}

kernel void ToulorgeDesmet2N84F_4a(global REAL *u1, global REAL *u2, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    REAL coef_A = (REAL)(-0.5515812888932000);

    REAL du[NUM_CONSERVED_VARS];
    compute_du_dt(time[1], ix, iy, iz, u1, du);

    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
        axpy_field_component(ix, iy, iz, i, coef_A, u2, (REAL)(DT) * du[i]);
    }
}
kernel void ToulorgeDesmet2N84F_4b(global REAL *u1, global REAL *u2, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    REAL coef_B = (REAL)(0.09911841297339970);
    REAL coef_C = (REAL)(0.3850364824285470);

    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
        axpy_field_component(ix, iy, iz, i, (REAL)(1), u1, coef_B * get_field_component(ix, iy, iz, i, u2));
    }

    if (ix == 0 && iy == 0 && iz == 0) {
        // time[0] is the time at the beginning of a time step
        // time[1] is the local time at one stage of a Runge-Kutta method
        time[1] = time[0] + coef_C * (REAL)(DT);
	}
}

kernel void ToulorgeDesmet2N84F_5a(global REAL *u1, global REAL *u2, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    REAL coef_A = (REAL)(-1.885790377558741);

    REAL du[NUM_CONSERVED_VARS];
    compute_du_dt(time[1], ix, iy, iz, u1, du);

    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
        axpy_field_component(ix, iy, iz, i, coef_A, u2, (REAL)(DT) * du[i]);
    }
}
kernel void ToulorgeDesmet2N84F_5b(global REAL *u1, global REAL *u2, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    REAL coef_B = (REAL)(0.7466920411064123);
    REAL coef_C = (REAL)(0.5040052477534100);

    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
        axpy_field_component(ix, iy, iz, i, (REAL)(1), u1, coef_B * get_field_component(ix, iy, iz, i, u2));
    }

    if (ix == 0 && iy == 0 && iz == 0) {
        // time[0] is the time at the beginning of a time step
        // time[1] is the local time at one stage of a Runge-Kutta method
        time[1] = time[0] + coef_C * (REAL)(DT);
	}
}

kernel void ToulorgeDesmet2N84F_6a(global REAL *u1, global REAL *u2, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    REAL coef_A = (REAL)(-5.701295742793264);

    REAL du[NUM_CONSERVED_VARS];
    compute_du_dt(time[1], ix, iy, iz, u1, du);

    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
        axpy_field_component(ix, iy, iz, i, coef_A, u2, (REAL)(DT) * du[i]);
    }
}
kernel void ToulorgeDesmet2N84F_6b(global REAL *u1, global REAL *u2, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    REAL coef_B = (REAL)(1.679584245618894);
    REAL coef_C = (REAL)(0.6578977561168540);

    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
        axpy_field_component(ix, iy, iz, i, (REAL)(1), u1, coef_B * get_field_component(ix, iy, iz, i, u2));
    }

    if (ix == 0 && iy == 0 && iz == 0) {
        // time[0] is the time at the beginning of a time step
        // time[1] is the local time at one stage of a Runge-Kutta method
        time[1] = time[0] + coef_C * (REAL)(DT);
	}
}

kernel void ToulorgeDesmet2N84F_7a(global REAL *u1, global REAL *u2, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    REAL coef_A = (REAL)(2.113903965664793);

    REAL du[NUM_CONSERVED_VARS];
    compute_du_dt(time[1], ix, iy, iz, u1, du);

    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
        axpy_field_component(ix, iy, iz, i, coef_A, u2, (REAL)(DT) * du[i]);
    }
}
kernel void ToulorgeDesmet2N84F_7b(global REAL *u1, global REAL *u2, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    REAL coef_B = (REAL)(0.2433728067008188);
    REAL coef_C = (REAL)(0.9484087623348481);

    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
        axpy_field_component(ix, iy, iz, i, (REAL)(1), u1, coef_B * get_field_component(ix, iy, iz, i, u2));
    }

    if (ix == 0 && iy == 0 && iz == 0) {
        // time[0] is the time at the beginning of a time step
        // time[1] is the local time at one stage of a Runge-Kutta method
        time[1] = time[0] + coef_C * (REAL)(DT);
	}
}

kernel void ToulorgeDesmet2N84F_8a(global REAL *u1, global REAL *u2, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    REAL coef_A = (REAL)(-0.5339578826675280);

    REAL du[NUM_CONSERVED_VARS];
    compute_du_dt(time[1], ix, iy, iz, u1, du);

    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
        axpy_field_component(ix, iy, iz, i, coef_A, u2, (REAL)(DT) * du[i]);
    }
}
kernel void ToulorgeDesmet2N84F_8b(global REAL *u1, global REAL *u2, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    REAL coef_B = (REAL)(0.1422730459001373);

    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
        axpy_field_component(ix, iy, iz, i, (REAL)(1), u1, coef_B * get_field_component(ix, iy, iz, i, u2));
    }

    if (ix == 0 && iy == 0 && iz == 0) {
        // time[0] is the time at the beginning of a time step
        // time[1] is the local time at one stage of a Runge-Kutta method
        time[0] = time[0] + (REAL)(DT);
        time[1] = time[0];
    }
}

//--------------------------------------------------------------------------------------------------
// Time integrator kernels using 3 fields
//--------------------------------------------------------------------------------------------------


/*
  KennedyCarpenterLewis2R54C_X

Perform one time step using the five stage, fourth order, explicit Runge-Kutta method
KennedyCarpenterLewis2R54C.
*/
kernel void KennedyCarpenterLewis2R54C_1a(global REAL *u1, global REAL *u2, global REAL *u3, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    REAL du[NUM_CONSERVED_VARS];
    compute_du_dt(time[1], ix, iy, iz, u1, du);

    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
        set_field_component(ix, iy, iz, i, u3, (REAL)(DT) * du[i]);
    }
}
kernel void KennedyCarpenterLewis2R54C_1b(global REAL *u1, global REAL *u2, global REAL *u3, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    REAL coef_A = (REAL)(0.22502245872571303);
    REAL coef_B = (REAL)(-0.17379315208537388);
    REAL coef_C = (REAL)(0.22502245872571303);

    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
        axpy_field_component(ix, iy, iz, i, (REAL)(1), u1, coef_A * get_field_component(ix, iy, iz, i, u3));
        aypz_field_component(ix, iy, iz, i, u2, 1, u1, coef_B * get_field_component(ix, iy, iz, i, u3));
    }

    if (ix == 0 && iy == 0 && iz == 0) {
        // time[0] is the time at the beginning of a time step
        // time[1] is the local time at one stage of a Runge-Kutta method
        time[1] = time[0] + coef_C * (REAL)(DT);
    }
}

kernel void KennedyCarpenterLewis2R54C_2a(global REAL *u1, global REAL *u2, global REAL *u3, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    REAL du[NUM_CONSERVED_VARS];
    compute_du_dt(time[1], ix, iy, iz, u1, du);

    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
        set_field_component(ix, iy, iz, i, u3, (REAL)(DT) * du[i]);
    }
}
kernel void KennedyCarpenterLewis2R54C_2b(global REAL *u1, global REAL *u2, global REAL *u3, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    REAL coef_A = (REAL)(0.5440433129514047);
    REAL coef_B = (REAL)(-0.1630884872250029);
    REAL coef_C = (REAL)(0.5952726195917439);

    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
        axpy_field_component(ix, iy, iz, i, (REAL)(1), u2, coef_A * get_field_component(ix, iy, iz, i, u3));
        aypz_field_component(ix, iy, iz, i, u1, 1, u2, coef_B * get_field_component(ix, iy, iz, i, u3));
    }

    if (ix == 0 && iy == 0 && iz == 0) {
        // time[0] is the time at the beginning of a time step
        // time[1] is the local time at one stage of a Runge-Kutta method
        time[1] = time[0] + coef_C * (REAL)(DT);
	}
}

kernel void KennedyCarpenterLewis2R54C_3a(global REAL *u1, global REAL *u2, global REAL *u3, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    REAL du[NUM_CONSERVED_VARS];
    compute_du_dt(time[1], ix, iy, iz, u2, du);

    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
        set_field_component(ix, iy, iz, i, u3, (REAL)(DT) * du[i]);
    }
}
kernel void KennedyCarpenterLewis2R54C_3b(global REAL *u1, global REAL *u2, global REAL *u3, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    REAL coef_A = (REAL)(0.14456824349399464);
    REAL coef_B = (REAL)(-0.5179208398863779);
    REAL coef_C = (REAL)(0.5767523758607357);

    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
        axpy_field_component(ix, iy, iz, i, (REAL)(1), u1, coef_A * get_field_component(ix, iy, iz, i, u3));
        aypz_field_component(ix, iy, iz, i, u2, 1, u1, coef_B * get_field_component(ix, iy, iz, i, u3));
    }
}

kernel void KennedyCarpenterLewis2R54C_4a(global REAL *u1, global REAL *u2, global REAL *u3, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    REAL du[NUM_CONSERVED_VARS];
    compute_du_dt(time[1], ix, iy, iz, u1, du);

    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
        set_field_component(ix, iy, iz, i, u3, (REAL)(DT) * du[i]);
    }
}
kernel void KennedyCarpenterLewis2R54C_4b(global REAL *u1, global REAL *u2, global REAL *u3, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    REAL coef_A = (REAL)(0.7866643421983568);
    REAL coef_B = (REAL)(-0.19416305717199442);
    REAL coef_C = (REAL)(0.8454958781727144);

    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
        axpy_field_component(ix, iy, iz, i, (REAL)(1), u2, coef_A * get_field_component(ix, iy, iz, i, u3));
        aypz_field_component(ix, iy, iz, i, u1, 1, u2, coef_B * get_field_component(ix, iy, iz, i, u3));
    }

    if (ix == 0 && iy == 0 && iz == 0) {
        // time[0] is the time at the beginning of a time step
        // time[1] is the local time at one stage of a Runge-Kutta method
        time[1] = time[0] + coef_C * (REAL)(DT);
	}
}

kernel void KennedyCarpenterLewis2R54C_5(global REAL *u1, global REAL *u2, global REAL *u3, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    REAL coef_B = (REAL)(0.34866717899927996);

    REAL du[NUM_CONSERVED_VARS];
    compute_du_dt(time[1], ix, iy, iz, u2, du);

    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
        axpy_field_component(ix, iy, iz, i, (REAL)(1), u1, coef_B * (REAL)(DT) * du[i]);
    }
}

/*
    CalvoFrancoRandez2R64_X

Perform one time step using the six stage, fourth order, explicit Runge-Kutta method
CalvoFrancoRandez2R64.
*/
kernel void CalvoFrancoRandez2R64_1a(global REAL *u1, global REAL *u2, global REAL *u3, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    REAL du[NUM_CONSERVED_VARS];
    compute_du_dt(time[1], ix, iy, iz, u1, du);

    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
        set_field_component(ix, iy, iz, i, u3, (REAL)(DT) * du[i]);
    }
}
kernel void CalvoFrancoRandez2R64_1b(global REAL *u1, global REAL *u2, global REAL *u3, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);


    REAL coef_A = (REAL)(0.17985400977138);
    REAL coef_B = (REAL)(0.10893125722541);
    REAL coef_C = (REAL)(0.28878526699679);

    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
        axpy_field_component(ix, iy, iz, i, (REAL)(1), u1, coef_B * get_field_component(ix, iy, iz, i, u3));
        aypz_field_component(ix, iy, iz, i, u2, 1, u1, coef_A * get_field_component(ix, iy, iz, i, u3));
    }

    if (ix == 0 && iy == 0 && iz == 0) {
        // time[0] is the time at the beginning of a time step
        // time[1] is the local time at one stage of a Runge-Kutta method
        time[1] = time[0] + coef_C * (REAL)(DT);
	}
}

kernel void CalvoFrancoRandez2R64_2a(global REAL *u1, global REAL *u2, global REAL *u3, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    REAL du[NUM_CONSERVED_VARS];
    compute_du_dt(time[1], ix, iy, iz, u2, du);

    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
        set_field_component(ix, iy, iz, i, u3, (REAL)(DT) * du[i]);
    }
}

kernel void CalvoFrancoRandez2R64_2b(global REAL *u1, global REAL *u2, global REAL *u3, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    REAL coef_A = (REAL)(0.14081893152111);
    REAL coef_B = (REAL)(0.13201701492152);
    REAL coef_C = (REAL)(0.38176720366804);

    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
        axpy_field_component(ix, iy, iz, i, (REAL)(1), u1, coef_B * get_field_component(ix, iy, iz, i, u3));
        aypz_field_component(ix, iy, iz, i, u2, 1, u1, coef_A * get_field_component(ix, iy, iz, i, u3));
    }

    if (ix == 0 && iy == 0 && iz == 0) {
        // time[0] is the time at the beginning of a time step
        // time[1] is the local time at one stage of a Runge-Kutta method
        time[1] = time[0] + coef_C * (REAL)(DT);
	}
}

kernel void CalvoFrancoRandez2R64_3a(global REAL *u1, global REAL *u2, global REAL *u3, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    REAL du[NUM_CONSERVED_VARS];
    compute_du_dt(time[1], ix, iy, iz, u2, du);

    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
        set_field_component(ix, iy, iz, i, u3, (REAL)(DT) * du[i]);
    };
}

kernel void CalvoFrancoRandez2R64_3b(global REAL *u1, global REAL *u2, global REAL *u3, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    REAL coef_A = (REAL)(0.08255631629428);
    REAL coef_B = (REAL)(0.38911623225517);
    REAL coef_C = (REAL)(0.71262082069639);

    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
        axpy_field_component(ix, iy, iz, i, (REAL)(1), u1, coef_B * get_field_component(ix, iy, iz, i, u3));
        aypz_field_component(ix, iy, iz, i, u2, 1, u1, coef_A * get_field_component(ix, iy, iz, i, u3));
    }

    if (ix == 0 && iy == 0 && iz == 0) {
        // time[0] is the time at the beginning of a time step
        // time[1] is the local time at one stage of a Runge-Kutta method
        time[1] = time[0] + coef_C * (REAL)(DT);
	}
}

kernel void CalvoFrancoRandez2R64_4a(global REAL *u1, global REAL *u2, global REAL *u3, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    REAL du[NUM_CONSERVED_VARS];
    compute_du_dt(time[1], ix, iy, iz, u2, du);

    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
        set_field_component(ix, iy, iz, i, u3, (REAL)(DT) * du[i]);
    }
}

kernel void CalvoFrancoRandez2R64_4b(global REAL *u1, global REAL *u2, global REAL *u3, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    REAL coef_A = (REAL)(0.65804425034331);
    REAL coef_B = (REAL)(-0.59203884581148);
    REAL coef_C = (REAL)(0.69606990893393);

    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
        axpy_field_component(ix, iy, iz, i, (REAL)(1), u1, coef_B * get_field_component(ix, iy, iz, i, u3));
        aypz_field_component(ix, iy, iz, i, u2, 1, u1, coef_A * get_field_component(ix, iy, iz, i, u3));
    }

    if (ix == 0 && iy == 0 && iz == 0) {
        // time[0] is the time at the beginning of a time step
        // time[1] is the local time at one stage of a Runge-Kutta method
        time[1] = time[0] + coef_C * (REAL)(DT);
	}
}

kernel void CalvoFrancoRandez2R64_5a(global REAL *u1, global REAL *u2, global REAL *u3, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    REAL du[NUM_CONSERVED_VARS];
    compute_du_dt(time[1], ix, iy, iz, u2, du);

    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
        set_field_component(ix, iy, iz, i, u3, (REAL)(DT) * du[i]);
    }
}

kernel void CalvoFrancoRandez2R64_5b(global REAL *u1, global REAL *u2, global REAL *u3, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    REAL coef_A = (REAL)(0.31862993413251);
    REAL coef_B = (REAL)(0.47385028714844);
    REAL coef_C = (REAL)(0.83050587987157);

    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
        axpy_field_component(ix, iy, iz, i, (REAL)(1), u1, coef_B * get_field_component(ix, iy, iz, i, u3));
        aypz_field_component(ix, iy, iz, i, u2, 1, u1, coef_A * get_field_component(ix, iy, iz, i, u3));
    }

    if (ix == 0 && iy == 0 && iz == 0) {
        // time[0] is the time at the beginning of a time step
        // time[1] is the local time at one stage of a Runge-Kutta method
        time[1] = time[0] + coef_C * (REAL)(DT);
	}
}

kernel void CalvoFrancoRandez2R64_6(global REAL *u1, global REAL *u2, global REAL *u3, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);


    REAL coef_B = (REAL)(0.48812405426094);

    REAL du[NUM_CONSERVED_VARS];
    compute_du_dt(time[1], ix, iy, iz, u2, du);

    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
        axpy_field_component(ix, iy, iz, i, (REAL)(1), u1, coef_B * (REAL)(DT) * du[i]);
    }
}

/*
  SSPRK33_X

Perform one time step using the three stage, third order, strong stability preserving explicit
Runge-Kutta method SSPRK33. Starting with the old value `b1`, the new value after one time step is
obtained using the temporary arrays `u2` and `u3`.
*/
kernel void SSPRK33_1(global REAL *u1, global REAL *u2, global REAL *u3, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    REAL du[NUM_CONSERVED_VARS];
    compute_du_dt(time[1], ix, iy, iz, u1, du);

    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
        aypz_field_component(ix, iy, iz, i, u2, 1, u1, (REAL)DT*du[i]);
    }
}

kernel void SSPRK33_2a(global REAL *u1, global REAL *u2, global REAL *u3, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    if (ix == 0 && iy == 0 && iz == 0) {
        // time[0] is the time at the beginning of a time step
        // time[1] is the local time at one stage of a Runge-Kutta method
        time[1] = time[0] + (REAL)(DT);
	  }
}
kernel void SSPRK33_2b(global REAL *u1, global REAL *u2, global REAL *u3, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    REAL du[NUM_CONSERVED_VARS];
    compute_du_dt(time[1], ix, iy, iz, u2, du);

    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
        aypz_field_component(ix, iy, iz, i, u3, (REAL)(0.75), u1, (REAL)(0.25) * (get_field_component(ix, iy, iz, i, u2) + (REAL)DT*du[i]));
    }
}

kernel void SSPRK33_3a(global REAL *u1, global REAL *u2, global REAL *u3, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    if (ix == 0 && iy == 0 && iz == 0) {
        // time[0] is the time at the beginning of a time step
        // time[1] is the local time at one stage of a Runge-Kutta method
        time[1] = time[0] + (REAL)(0.5*DT);
	  }
}
kernel void SSPRK33_3b(global REAL *u1, global REAL *u2, global REAL *u3, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    REAL du[NUM_CONSERVED_VARS];
    compute_du_dt(time[1], ix, iy, iz, u3, du);

    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
        aypz_field_component(ix, iy, iz, i, u1, (REAL)(1.0/3.0), u1, (REAL)(2.0/3.0) * (get_field_component(ix, iy, iz, i, u3) + (REAL)DT*du[i]));
    }
}

/*
  SSPRK104_X

Perform one time step using the ten stage, fourth order, strong stability preserving explicit
Runge-Kutta method SSPRK104.

Note: SSPRK104_06 is just a linear combination step and does not require the computation of a time
derivative.
*/
kernel void SSPRK104_01(global REAL *u1, global REAL *u2, global REAL *u3, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    REAL du[NUM_CONSERVED_VARS];
    compute_du_dt(time[1], ix, iy, iz, u1, du);

    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
        aypz_field_component(ix, iy, iz, i, u2, 1, u1, (REAL)(DT/6.0) * du[i]);
    }
}

kernel void SSPRK104_02a(global REAL *u1, global REAL *u2, global REAL *u3, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    if (ix == 0 && iy == 0 && iz == 0) {
        // time[0] is the time at the beginning of a time step
        // time[1] is the local time at one stage of a Runge-Kutta method
        time[1] = time[0] + (REAL)(DT/6.0);
	  }
}
kernel void SSPRK104_02b(global REAL *u1, global REAL *u2, global REAL *u3, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    REAL du[NUM_CONSERVED_VARS];
    compute_du_dt(time[1], ix, iy, iz, u2, du);

    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
        aypz_field_component(ix, iy, iz, i, u3, 1, u2, (REAL)(DT/6.0) * du[i]);
    }
}

kernel void SSPRK104_03a(global REAL *u1, global REAL *u2, global REAL *u3, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    if (ix == 0 && iy == 0 && iz == 0) {
        // time[0] is the time at the beginning of a time step
        // time[1] is the local time at one stage of a Runge-Kutta method
        time[1] = time[0] + (REAL)(DT/3.0);
	  }
}
kernel void SSPRK104_03b(global REAL *u1, global REAL *u2, global REAL *u3, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    REAL du[NUM_CONSERVED_VARS];
    compute_du_dt(time[1], ix, iy, iz, u3, du);

    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
        aypz_field_component(ix, iy, iz, i, u2, 1, u3, (REAL)(DT/6.0) * du[i]);
    }
}

kernel void SSPRK104_04a(global REAL *u1, global REAL *u2, global REAL *u3, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);


    if (ix == 0 && iy == 0 && iz == 0) {
        // time[0] is the time at the beginning of a time step
        // time[1] is the local time at one stage of a Runge-Kutta method
        time[1] = time[0] + (REAL)(0.5*DT);
	  }
}
kernel void SSPRK104_04b(global REAL *u1, global REAL *u2, global REAL *u3, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    REAL du[NUM_CONSERVED_VARS];
    compute_du_dt(time[1], ix, iy, iz, u2, du);

    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
        aypz_field_component(ix, iy, iz, i, u3, 1, u2, (REAL)(DT/6.0) * du[i]);
    }
}

kernel void SSPRK104_05a(global REAL *u1, global REAL *u2, global REAL *u3, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    if (ix == 0 && iy == 0 && iz == 0) {
        // time[0] is the time at the beginning of a time step
        // time[1] is the local time at one stage of a Runge-Kutta method
        time[1] = time[0] + (REAL)(DT*2.0/3.0);
	  }
}
kernel void SSPRK104_05b(global REAL *u1, global REAL *u2, global REAL *u3, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    REAL du[NUM_CONSERVED_VARS];
    compute_du_dt(time[1], ix, iy, iz, u3, du);

    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
        aypz_field_component(ix, iy, iz, i, u2, 1, u3, (REAL)(DT/6.0) * du[i]);
    }
}

kernel void SSPRK104_06(global REAL *u1, global REAL *u2, global REAL *u3, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
        aypz_field_component(ix, iy, iz, i, u3, (REAL)(0.04), u1, (REAL)(0.36) * get_field_component(ix, iy, iz, i, u2));
        aypz_field_component(ix, iy, iz, i, u2, (REAL)(0.60), u1, (REAL)(0.40) * get_field_component(ix, iy, iz, i, u2));
    }

    if (ix == 0 && iy == 0 && iz == 0) {
        // time[0] is the time at the beginning of a time step
        // time[1] is the local time at one stage of a Runge-Kutta method
        time[1] = time[0] + (REAL)(DT/3.0);
	  }
}

kernel void SSPRK104_07(global REAL *u1, global REAL *u2, global REAL *u3, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    REAL du[NUM_CONSERVED_VARS];
    compute_du_dt(time[1], ix, iy, iz, u2, du);

    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
        aypz_field_component(ix, iy, iz, i, u1, 1, u2, (REAL)(DT/6.0) * du[i]);
    }
}

kernel void SSPRK104_08a(global REAL *u1, global REAL *u2, global REAL *u3, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    if (ix == 0 && iy == 0 && iz == 0) {
        // time[0] is the time at the beginning of a time step
        // time[1] is the local time at one stage of a Runge-Kutta method
        time[1] = time[0] + (REAL)(0.5*DT);
  	}
}
kernel void SSPRK104_08b(global REAL *u1, global REAL *u2, global REAL *u3, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    REAL du[NUM_CONSERVED_VARS];
    compute_du_dt(time[1], ix, iy, iz, u1, du);

    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
        aypz_field_component(ix, iy, iz, i, u2, 1, u1, (REAL)(DT/6.0) * du[i]);
    }
}

kernel void SSPRK104_09a(global REAL *u1, global REAL *u2, global REAL *u3, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    if (ix == 0 && iy == 0 && iz == 0) {
        // time[0] is the time at the beginning of a time step
        // time[1] is the local time at one stage of a Runge-Kutta method
        time[1] = time[0] + (REAL)(DT*2.0/3.0);
  	}
}
kernel void SSPRK104_09b(global REAL *u1, global REAL *u2, global REAL *u3, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    REAL du[NUM_CONSERVED_VARS];
    compute_du_dt(time[1], ix, iy, iz, u2, du);

    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
        aypz_field_component(ix, iy, iz, i, u1, 1, u2, (REAL)(DT/6.0) * du[i]);
    }
}

kernel void SSPRK104_10a(global REAL *u1, global REAL *u2, global REAL *u3, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    if (ix == 0 && iy == 0 && iz == 0) {
        // time[0] is the time at the beginning of a time step
        // time[1] is the local time at one stage of a Runge-Kutta method
        time[1] = time[0] + (REAL)(DT*5.0/6.0);
  	}
}
kernel void SSPRK104_10b(global REAL *u1, global REAL *u2, global REAL *u3, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    REAL du[NUM_CONSERVED_VARS];
    compute_du_dt(time[1], ix, iy, iz, u1, du);

    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
        aypz_field_component(ix, iy, iz, i, u2, 1, u1, (REAL)(DT/6.0) * du[i]);
    }
}

kernel void SSPRK104_11a(global REAL *u1, global REAL *u2, global REAL *u3, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    if (ix == 0 && iy == 0 && iz == 0) {
        // time[0] is the time at the beginning of a time step
        // time[1] is the local time at one stage of a Runge-Kutta method
        time[0] = time[0] + (REAL)(DT);
        time[1] = time[0];
	  }
}
kernel void SSPRK104_11b(global REAL *u1, global REAL *u2, global REAL *u3, global REAL *time) {

    uint ix = get_global_id(0);
    uint iy = get_global_id(1);
    uint iz = get_global_id(2);

    REAL du[NUM_CONSERVED_VARS];
    compute_du_dt(time[1], ix, iy, iz, u2, du);

    for (uint i = 0; i < NUM_CONSERVED_VARS; ++i) {
        aypz_field_component(ix, iy, iz, i, u1, 1, u3, (REAL)(0.6) * (get_field_component(ix, iy, iz, i, u2) + (REAL)(DT/6.0) * du[i]));
    }
}
