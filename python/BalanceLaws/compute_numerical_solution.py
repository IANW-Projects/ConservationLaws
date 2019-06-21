def compute_numerical_solution(field_u1, field_u2, c_time = 'none'):
    global I_mesh, I_TI, I_BalanceLaws, I_Tech, I_RunOps, I_Results
    
    
    if I_Tech['REAL'] == 'float':
        current_time = np.zeros(2, dtype=np.float32)
        components = np.zeros(2, dtype=np.float32)
        norm_output = np.zeros((I_Tech['num_groups'], 1), dtype=np.float32)
        Lerror = np.zeros((I_TI['num_steps']+1, I_BalanceLaws['NUM_CONSERVED_VARS']), dtype=np.float32)
    else:
        current_time = np.zeros(2, dtype=np.float64)
        components = np.zeros(2, dtype=np.float64)
        norm_output = np.zeros((I_Tech['num_groups'], 1), dtype=np.float64)
        Lerror = np.zeros((I_TI['num_steps']+1, I_BalanceLaws['NUM_CONSERVED_VARS']), dtype=np.float64)
    
    if c_time != 'none'
        current_time = c_time

    if I_TI['time_integrator'] == 'CarpenterKennedy2N54':
        RK_Step_a = ['CarpenterKennedy2N54_1a', 'CarpenterKennedy2N54_2a', \
                     'CarpenterKennedy2N54_3a', 'CarpenterKennedy2N54_4a', \
                     'CarpenterKennedy2N54_5a']

        RK_Step_b = ['CarpenterKennedy2N54_1b', 'CarpenterKennedy2N54_2b', \
                     'CarpenterKennedy2N54_3b', 'CarpenterKennedy2N54_4b', \
                     'CarpenterKennedy2N54_5b']
        RK_Step = []
        for i in range(len(RK_Step_a)):
            RK_Step.append('calc_auxiliary_vars_1_2_args')
            RK_Step.append(RK_Step_a[i])
            RK_Step.append(RK_Step_b[i])

        time_integrator_num_fields = 2;

    if I_TI['time_integrator'] == 'SSPRK33':
        RK_Step = []
        RK_Step.append('calc_auxiliary_vars_1_3_args')
        RK_Step.append('SSPRK33_1')
        RK_Step.append('SSPRK33_2a')
        RK_Step.append('calc_auxiliary_vars_2_3_args')
        RK_Step.append('SSPRK33_2b')
        RK_Step.append('SSPRK33_3a')
        RK_Step.append('calc_auxiliary_vars_3_3_args')
        RK_Step.append('SSPRK33_3b')
        RK_Step.append('calc_time_3_args')
        time_integrator_num_fields = 3

    num_steps_run = I_TI['num_steps']
    kernel_runtime = 0
    RK_block_size = 15000

    if time_integrator_num_fields == 2:
        if I_RunOps['save_integrals_over_time']:
            for step in range(I_TI['num_steps']):
                for comp in range(I_BalanceLaws['NUM_CONSERVED_VARS']):
                    components[0] = comp
                    cl.run_kernel('analytical_u', I_BalanceLaws['g_range'], I_BalanceLaws['l_range'], field_u2, current_time)
                    norm_output[:] = 0
                    if I_RunOps['norm'] == 'L2':
                            cl.run_kernel('norm2_diff'. I_Tech['g_range'], I_Tech['l_range'], field_u1, field_u2, norm_output, components)
                            Lerror[comp, step] = sqrt(np.sum(norm_output))
                    elif I_RunOps['norm'] == 'LInf':
                            cl.run_kernel('norm_infty_diff', I_Tech['g_range'], I_Tech['l_range'], field_u1, field_u2, norm_output, components)
                            Lerror[comp, step] = np.max(norm_output)
                cl.run_kernel(RK_Step, I_BalanceLaws['g_range'], I_BalanceLaws['l_range'], field_u1, field_u2, current_time)
        else:
            while num_steps_run > RK_block_size:
                kernel_list = RK_Step * RK_block_size
                t = cl.run_kernel(kernel_list, I_BalanceLaws['g_range'], I_BalanceLaws['l_range'], field_u1, field_u2, current_time)
                kernel_runtime = kernel_runtime + t
                num_steps_run = num_steps_run - RK_block_size
            if num_steps_run > 0:
                kernel_list  = num_steps_run* RK_Step
                t = cl.run_kernel(kernel_list,  I_BalanceLaws['g_range'], I_BalanceLaws['l_range'], field_u1, field_u2, current_time)
                kernel_runtime = kernel_runtime + t


    if time_integrator_num_fields == 3:
        num_nodes = I_Mesh['NODES_X']*I_Mesh['NODES_Y']*I_Mesh['NODES_Z']
        if I_Tech['REAL'] == 'float':
            field_u3 = np.zeros((num_nodes*I_BalanceLaws['NUM_TOTAL_VARS'], 1), dtype = np.float32)
        else:
            field_u3 = np.zeros((num_nodes*I_BalanceLaws['NUM_TOTAL_VARS'], 1), dtype = np.float64)

        if I_RunOps['save_integrals_over_time']:
            for step in range(TI['num_steps']):
                for comp in range(I_BalanceLaws['NUM_CONSERVED_VARS']):
                    components[0] = comp
                    cl.run_kernel('analytical_u', I_BalanceLaws['g_range'], I_BalanceLaws['l_range'], field_u2, current_time)
                    norm_output[:] = 0
                    if I_RunOps['norm'] == 'L2':
                            cl.run_kernel('norm2_diff'. I_Tech['g_range'], I_Tech['l_range'], field_u1, field_u2, field_u3,norm_output, components)
                            Lerror[comp, step] = sqrt(np.sum(norm_output))
                    elif I_RunOps['norm'] == 'LInf':
                            cl.run_kernel('norm_infty_diff', I_Tech['g_range'], I_Tech['l_range'], field_u1, field_u2, field_u3,norm_output, components)
                            Lerror[comp, step] = np.max(norm_output)
                cl.run_kernel(RK_Step, I_BalanceLaws['g_range'], I_BalanceLaws['l_range'], field_u1, field_u2, field_u3, current_time)

        else:
            while num_steps_run > RK_block_size:
                kernel_list = RK_Step * RK_block_size
                t = cl.run_kernel(kernel_list, I_BalanceLaws['g_range'], I_BalanceLaws['l_range'], field_u1, field_u2, field_u3, current_time)
                kernel_runtime = kernel_runtime + t
                num_steps_run = num_steps_run - RK_block_size
            if num_steps_run > 0:
                kernel_list  = num_steps_run* RK_Step
                t = cl.run_kernel(kernel_list,  I_BalanceLaws['g_range'], I_BalanceLaws['l_range'], field_u1, field_u2, field_u3, current_time)
                kernel_runtime = kernel_runtime + t


    if I_RunOps['save_fields']:
        I_Results['field_u'] = field_u1

    print(f"current_time[0] = {current_time[0]} current_time[1] = {current_time[1]}")
    current_time[1] = I_TI['final_time']
    cl.run_kernel('analytical_u', I_BalanceLaws['g_range'], I_BalanceLaws['l_range'], field_u2, current_time)
    abs_err = np.zeros((I_BalanceLaws['NUM_CONSERVED_VARS'],1), dtype = np.float64)
    rel_err = np.zeros(( I_BalanceLaws['NUM_CONSERVED_VARS'],1),dtype=np.float64)

    for comp in range(I_BalanceLaws['NUM_CONSERVED_VARS']):
        components[0] = comp
        norm_output[:]= 0
        if I_RunOps['norm'] ==  'L2':
#            cl.run_kernel('norm2_diff', I_Tech['g_range'], I_Tech['l_range'],field_u1, field_u2, norm_output, components)
            abs_err[comp] = math.sqrt(np.sum(norm_output))
            norm_output[:] = 0
  #          cl.run_kernel('norm2', I_Tech['g_range'], I_Tech['l_range'], field_u2, norm_output, components)
            rel_err[comp] = abs_err[comp] / math.sqrt(np.sum(norm_output))
        elif I_RunOps['norm'] == 'LInf':
#            cl.run_kernel('norm_infty_diff', I_Tech['g_range'], I_Tech['l_range'],field_u1, field_u2, norm_output, components)
            abs_err[comp] = np.max(norm_output)
            norm_output[:] = 0
 #           cl.run_kernel('norm_infty', I_Tech['g_range'], I_Tech['l_range'], field_u2, norm_output, components)
            rel_err[comp] = abs_err[comp] / np.max(norm_output)

    I_Results['abs_err'] = abs_err
    I_Results['rel_err'] = rel_err

    if I_RunOps['save_integrals_over_time']:
        Lerror[I_TI['num_steps'],:] = I_Results['abs_err']
        I_Results['Lerror_over_time'] = Lerror
        I_Results['time'] = np.linspace(0, I_TI['final_time'], len(Lerror)*2)
    return current_time
