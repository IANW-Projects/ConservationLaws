%This project is licensed under the terms of the Creative Commons CC BY-NC-ND 4.0 license.

function [] = prepare_vars()

    %prepare_vars() initializes containers (maps) in which all variables (label, value) are
    %stored. Variables are grouped together by purpose:
    %I_Mesh: Contains all variables concerning the discrete mesh, e.g.
    %stepsize, box size
    %I_TI: Contains all variables concerning time integration, e.g. time
    %integrator, timestep
    %I_BalanceLaws: Contains all variables relevant for the BalanceLaws equations, e.g.
    %discretization form, switch for additional artificial dissipation
    %I_Tech: Contains all variables of technical nature, e.g. the device,
    %computation precision
    %I_RunOps: Specifies the parameters of a computation, e.g. which variables
    %will be saved, the testcase
    %Variables in capital letters are program specific settings which are set using
    %OpenCL compiler defines. If the value is a string in captial letters
    %beginning with USE it acts as a switch, e.g. to enable artificial dissipation
    %or to switch between different discretizations. In case the key is written
    %in capital letters the program define will be set to the corresponding
    %value of the key.

    global I_Mesh I_TI I_BalanceLaws I_Tech I_RunOps I_Results

    keySet = {'NODES_X', 'NODES_Y', 'NODES_Z', 'DX', 'DY', 'DZ', 'XMIN','YMIN', 'ZMIN', 'XMAX', 'YMAX','ZMAX'};
    valueSet = {uint16(0) uint16(0) uint16(0) 0 0 0 0 0 0 0 0 0};
    I_Mesh = containers.Map(keySet, valueSet, 'UniformValues',false);

    keySet = {'cfl', 'final_time', 'time_integrator', 'DT', 'num_steps'};
    valueSet = {0 0 '' 0 0};
    I_TI = containers.Map(keySet, valueSet,'UniformValues',false);

    keySet = {'hall_term', 'g_range', 'l_range','NUM_CONSERVED_VARS', 'NUM_AUXILIARY_VARS', 'NUM_TOTAL_VARS'};
    valueSet = {'' 0 0 0 0 0};
    I_BalanceLaws = containers.Map(keySet, valueSet,'UniformValues',false);

    keySet = {'device', 'REAL', 'optimizations', 'NUM_NODES_PAD', 'num_groups','W_SIZE', 'g_range', 'l_range', 'memory_layout'};
    valueSet = {0 '' '' 0 0 0 0 0 'USE_STRUCTURE_OF_ARRAYS'};
    I_Tech = containers.Map(keySet, valueSet,'UniformValues',false);

    keySet = {'order', 'operator_form','conservation_laws', 'testcase', 'periodic', 'plot_numerical_solution', 'save_fields', 'save_integrals_over_time', 'norm'};
    valueSet = {0 'classical' '' '' '' '' false false 'L2'};
    I_RunOps = containers.Map(keySet, valueSet,'UniformValues',false);

    keySet = {'abs_err', 'rel_err', 'field_u', 'runtime', 'kernel_runtime' 'error_over_time' 'time'};
    valueSet = {0 0 0 0 0 0 0};
    I_Results = containers.Map(keySet, valueSet,'UniformValues',false);


end
