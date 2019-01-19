//This project is licensed under the terms of the Creative Commons CC BY-NC-ND 4.0 license.

// Contains kernel for norm and derivative operators
//--------------------------------------------------------------------------------------------------
// Norm
//--------------------------------------------------------------------------------------------------

/*
Computes the squared L^2 norm of the difference of the fields `d_field_1`, `d_field_2` for the field component `comp` and
stores the result in `output`.
*/
kernel void norm2_diff(global REAL * u1, global REAL *u2, global REAL *output, global REAL *comp) {

    const int gid = get_global_id(0);
    const int lid = get_local_id(0);
    const int group_size = get_local_size(0);
    const int wgid = get_group_id(0);

    local REAL partial_n[(uint)W_SIZE];

    uint4 s_idx = calc_sub_idx(gid);

    int bound_x = get_bound_x(s_idx.x, NUM_BOUNDS);
    int bound_y = get_bound_y(s_idx.y, NUM_BOUNDS);
    int bound_z = get_bound_z(s_idx.z, NUM_BOUNDS);

    REAL fac = ((REAL)DX / M_INV[NUM_BOUNDS + bound_x])
             * ((REAL)DY / M_INV[NUM_BOUNDS + bound_y])
             * ((REAL)DZ / M_INV[NUM_BOUNDS + bound_z]);

    REAL diff = get_field_component(s_idx.x, s_idx.y, s_idx.z, (uint)comp[0], u1) - get_field_component(s_idx.x, s_idx.y, s_idx.z, (uint)comp[0], u2);
    partial_n[lid] = fac * diff * diff;

    barrier(CLK_LOCAL_MEM_FENCE);

    for (uint i = group_size/2; i > 0; i /= 2) {
        if(lid < i) {
        partial_n[lid] += partial_n[lid + i];
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    if (lid == 0) {
        #ifdef REDUCE
        atomic_add_global(&(output[0]), partial_n[0]);
        #else
        output[wgid] = partial_n[0];
        #endif
    }
}


/*
Computes the squared L^2 norm for the field component `comp` of the field `d_field` and stores the result in `output`.
*/
kernel void norm2(global REAL *u, global REAL *output, global REAL* comp) {

  const int gid = get_global_id(0);
  const int lid = get_local_id(0);
  const int group_size = get_local_size(0);
  const int wgid = get_group_id(0);

  local REAL partial_n[(uint)W_SIZE];

  uint4 s_idx = calc_sub_idx(gid);

  int bound_x = get_bound_x(s_idx.x, NUM_BOUNDS);
  int bound_y = get_bound_y(s_idx.y, NUM_BOUNDS);
  int bound_z = get_bound_z(s_idx.z, NUM_BOUNDS);

  REAL fac = ((REAL)DX / M_INV[NUM_BOUNDS + bound_x])
           * ((REAL)DY / M_INV[NUM_BOUNDS + bound_y])
           * ((REAL)DZ / M_INV[NUM_BOUNDS + bound_z]);


  partial_n[lid] = fac * get_field_component(s_idx.x, s_idx.y, s_idx.z, (uint)comp[0], u)*get_field_component(s_idx.x, s_idx.y, s_idx.z, (uint)comp[0], u);

  barrier(CLK_LOCAL_MEM_FENCE);

  for (uint i = group_size/2; i > 0; i /= 2) {
    if(lid < i) {
      partial_n[lid] += partial_n[lid + i];
    }
    barrier(CLK_LOCAL_MEM_FENCE);
  }

  if (lid == 0) {
    #ifdef REDUCE
      atomic_add_global(&(output[0]), partial_n[0]);
    #else
      output[wgid] = partial_n[0];
    #endif
  }
}

/*
Computes the L^Infinity norm of the difference of the fields `d_field_1`, `d_field_2` for the field component `comp` and
stores the result in `output`.
*/

kernel void norm_infty_diff(global REAL * u1, global REAL *u2, global REAL *output, global REAL* comp) {

  const int gid = get_global_id(0);
  const int lid = get_local_id(0);
  const int group_size = get_local_size(0);
  const int wgid = get_group_id(0);

  local REAL local_max[(uint)W_SIZE];

  uint4 s_idx = calc_sub_idx(gid);

  REAL pmax = fabs(get_field_component(s_idx.x, s_idx.y, s_idx.z, (uint)comp[0], u1) - get_field_component(s_idx.x, s_idx.y, s_idx.z, (uint)comp[0], u2));

  if (lid == 0) {
    local_max[0] = (REAL)(0);
  }

  for (uint i = 0; i < group_size; i++) {
    barrier(CLK_LOCAL_MEM_FENCE);
    if(lid == i) {
      local_max[0] = max(local_max[0], pmax);
    }
  }

  if (lid == 0) {
    output[wgid] = local_max[0];
  }
}

/*
Computes the L^Infinity norm for the field component `comp` of the field `d_field` and stores the result in `output`.
*/

kernel void norm_infty(global REAL *u, global REAL *output, global REAL* comp) {

  const int gid = get_global_id(0);
  const int lid = get_local_id(0);
  const int group_size = get_local_size(0);
  const int wgid = get_group_id(0);

  local REAL local_max[(uint)W_SIZE];

  uint4 s_idx = calc_sub_idx(gid);

  REAL pmax = fabs(get_field_component(s_idx.x, s_idx.y, s_idx.z, (uint)comp[0], u));

  if (lid == 0) {
    local_max[0] = (REAL)(0);
  }

  for (uint i = 0; i < group_size; i++) {
    barrier(CLK_LOCAL_MEM_FENCE);
    if(lid == i) {
      local_max[0] = max(local_max[0], pmax);
    }
  }

  if (lid == 0) {
    output[wgid] = local_max[0];
  }
}

kernel void analytical_u(global REAL *u, global REAL *time) {

   uint ix = get_global_id(0);
   uint iy = get_global_id(1);
   uint iz = get_global_id(2);

   analytical_solution(ix, iy, iz, u, time[0]);
}
