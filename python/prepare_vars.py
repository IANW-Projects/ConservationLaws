# from +Balance_Laws/prepare_vars.m


def prepare_vars():
    global I_Mesh, I_TI, I_BalanceLaws, I_Tech, I_RunOps, I_Results
    I_Mesh = {'NODES_X':0, 'NODES_Y':0, 'NODES_Z':0, 'DX':0, 'DY':0, 'DZ':0, 'XMIN':0, 'XMAX':0, 'YMAX':0, 'YMIN':0, 'ZMIN':0, 'ZMAX':0}

    I_TI = {'cfl':0, 'final_time':0, 'time_integrator':'\0', 'DT':0, 'num_steps':0}

    I_Tech = {'device':0, 'REAL':'\0','optimizations':'\0','NUM_NODES_PAD':0, 'num_groups':0,'W_SIZE':0, 'g_range':0, 'l_range':0}

    I_BalanceLaws = {'hall_term':'\0', 'g_range':0, 'l_range':0, 'NUM_CONSERVED_VARS':0, 'NUM_AUXILIARY_VARS':0, 'NUM_TOTAL_VARS':0}
    
    I_RunOps = {'order':0, 'operator_form':'classical', 'conservation_laws':'\0', 'testcase':'\0', 'periodicx':'\0', 'periodicy':'\0', 'periodicz':'\0', 'plot_numerical_solution':'\0',  'save_fields':False, 'save_integrals_over_time':False, 'norm':'L2'}

    I_Results = {'abs_err':0, 'rel_err':0, 'field_u':0,  'runtime':0, 'kernel_runtime':0, 'error_over_time':0, 'time':0}

