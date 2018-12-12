import cl

def initialize():
    global I_Mesh, I_TI, I_BalanceLaws, I_Tech, I_RunOps
    
    cl.create_ctx()
        
    # a python float is a 64bit float => float means double
    if I_RunOps['periodic'] == 'USE_PERIODIC':
        DX = float(I_Mesh['XMAX'] - I_Mesh['XMIN']) / float(I_Mesh['NODES_X'])
        DY = float(I_Mesh['YMAX'] - I_Mesh['YMIN']) / float(I_Mesh['NODES_Y'])
        DZ = float(I_Mesh['ZMAX'] - I_Mesh['ZMIN']) / float(I_Mesh['NODES_Z'])
    else:
        DX = float(I_Mesh['XMAX'] - I_Mesh['XMIN']) / (float(I_Mesh['NODES_X'])-1)
        DY = float(I_Mesh['YMAX'] - I_Mesh['YMIN']) / (float(I_Mesh['NODES_Y'])-1)
        DZ = float(I_Mesh['ZMAX'] - I_Mesh['ZMIN']) / (float(I_Mesh['NODES_Z'])-1)

    I_Mesh['DX'] = DX; I_Mesh['DY'] = DY; I_Mesh['DZ'] = DZ

    I_BalanceLaws['g_range'] = np.array([I_Mesh['NODES_X'], I_Mesh['NODES_Y'], I_Mesh['NODES_Z']], dtype=np.uint32)
    I_BalanceLaws['l_range'] = np.array([0], dtype=np.uint32)

    num_nodes = I_Mesh['NODES_X']*I_Mesh['NODES_Y']*I_Mesh['NODES_Z']

    if (num_nodes > cl.group_size):
        group_size = cl.group_size
    else:
        group_size = 2**floor(log(num_nodes) / log(2))

    num_nodes_pad = math.ceil(num_nodes/group_size)*group_size
    num_groups = math.ceil(num_nodes_pad/group_size)

    I_Tech['num_nodes_pad'] = int(num_nodes_pad)
    I_Tech['num_groups'] = int(num_groups)
    I_Tech['W_SIZE'] = int(group_size)

    I_Tech['g_range'] = [I_Tech['num_nodes_pad'], 1,1]  
    I_Tech['l_range' ] = [I_Tech['W_SIZE'], 1,1]


    if I_Tech['REAL'] == 'float':
        field_u1 = np.zeros((num_nodes_pad*I_BalanceLaws['NUM_TOTAL_VARS'], 1), dtype=np.float32)
        field_u2 = np.zeros((num_nodes_pad*I_BalanceLaws['NUM_TOTAL_VARS'], 1), dtype=np.float32)
    else:
        field_u1 = np.zeros((num_nodes_pad*I_BalanceLaws['NUM_TOTAL_VARS'], 1), dtype = np.float64)
        field_u2 = np.zeros((num_nodes_pad*I_BalanceLaws['NUM_TOTAL_VARS'], 1), dtype = np.float64)


    x = np.linspace(I_Mesh['XMIN'], I_Mesh['XMAX'], I_Mesh['NODES_X']);
    y = np.linspace(I_Mesh['YMIN'], I_Mesh['YMAX'], I_Mesh['NODES_Y']);
    z = np.linspace(I_Mesh['ZMIN'], I_Mesh['ZMAX'], I_Mesh['NODES_Z']);

    kernel_path_list =  []

    if I_RunOps['order'] == 2:
        if I_RunOps['operator_form'] == 'extended':
            kernel_path_list.append('../include/2ndOrderExtended.h')
        else: # classical operators
            kernel_path_list.append('../include/2ndOrder.h')

    elif I_RunOps['order'] == 4:
        if I_RunOps['operator_form'] ==  'extended':
            kernel_path_list.append('../include/4thOrderExtended.h')
        else: # classical operators
            kernel_path_list.append('../include/4thOrder.h')

    elif I_RunOps['order'] == 6:
        if I_RunOps['operator_form'] == 'extended':
            kernel_path_list.append('../include/6thOrderExtended.h')
        else: # classical operators
            kernel_path_list.append('../include/6thOrder.h')

    else:
            print('Specify order \n')

    kernel_path_list.append('../include_defines/' +  I_RunOps['conservation_laws'] + '_defines.h')
    kernel_path_list.append('../include/utils.h')
    kernel_path_list.append('../kernel/SBP_operator.cl')
    kernel_path_list.append('../include_testcases/' +  I_RunOps['conservation_laws'] + '_' +  I_RunOps['testcase']+ '.h')
    kernel_path_list.append('../include_physics/' + I_RunOps['conservation_laws'] + '.h')
    kernel_path_list.append('../kernel/kernel_init.cl')

    settings_tech = generate_settings(I_Tech, ['REAL', 'REAL4', 'optimizations'])

    settings_mesh = generate_settings(I_Mesh,  ['DX', 'DY', 'DZ','NODES_X', 'NODES_Y', 'NODES_Z', 'XMIN', 'XMAX', 'YMIN', 'YMAX','ZMIN', 'ZMAX'])
    settings_runops = generate_settings(I_RunOps, ['periodic'])

    settings = settings_tech + settings_mesh + settings_runops

    cl.compile_kernels(kernel_path_list, settings)

    cl.run_kernel( 'init', I_BalanceLaws['g_range'], I_BalanceLaws['l_range'], field_u1)


    kernel_path_list =  []

    if I_RunOps['order'] == 2:
        if I_RunOps['operator_form'] == 'extended':
            kernel_path_list.append('../include/2ndOrderExtended.h')
        else: # classical operators
            kernel_path_list.append('../include/2ndOrder.h')

    elif I_RunOps['order'] == 4:
        if I_RunOps['operator_form'] ==  'extended':
            kernel_path_list.append('../include/4thOrderExtended.h')
        else: # classical operators
            kernel_path_list.append('../include/4thOrder.h')

    elif I_RunOps['order'] == 6:
        if I_RunOps['operator_form'] == 'extended':
            kernel_path_list.append('../include/6thOrderExtended.h')
        else: # classical operators
            kernel_path_list.append('../include/6thOrder.h')

    else:
            print('Specify order \n')

    kernel_path_list.append('../include_defines/' +  I_RunOps['conservation_laws'] + '_defines.h')
    kernel_path_list.append('../include/utils.h')
    kernel_path_list.append('../kernel/SBP_operator.cl')
    kernel_path_list.append('../include_testcases/' +  I_RunOps['conservation_laws'] + '_' +  I_RunOps['testcase']+ '.h')
    kernel_path_list.append('../include_physics/' + I_RunOps['conservation_laws'] + '.h')
    kernel_path_list.append('../kernel/kernel_time_integrator.cl')
    kernel_path_list.append('../kernel/kernel_norm.cl')


    settings_tech = generate_settings(I_Tech, ['REAL', 'REAL4', 'W_SIZE', 'optimizations'])

    settings_mesh = generate_settings(I_Mesh,  ['DX', 'DY', 'DZ','NODES_X', 'NODES_Y', 'NODES_Z', 'XMIN', 'XMAX', 'YMIN', 'YMAX','ZMIN', 'ZMAX'])
    settings_time_integration = generate_settings(I_TI, ['DT'])
    settings_runops = generate_settings(I_RunOps, ['periodic'])

    settings = settings_tech + settings_mesh + settings_time_integration +settings_runops

    cl.compile_kernels(kernel_path_list, settings)

    return field_u1, field_u2



