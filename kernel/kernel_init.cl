//This project is licensed under the terms of the Creative Commons CC BY-NC-ND 4.0 license.

kernel void init(global REAL *u) {

   uint ix = get_global_id(0);
   uint iy = get_global_id(1);
   uint iz = get_global_id(2);

   init_fields(ix, iy, iz, u);

}

