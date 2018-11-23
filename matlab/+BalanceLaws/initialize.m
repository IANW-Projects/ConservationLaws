%This project is licensed under the terms of the Creative Commons CC BY-NC-ND 4.0 license.

function [field_u1, field_u2] = initialize()

    %Initialize the magnetic field, the velocity field and the density field
    %according to the specified testcase. Also calculates and sets additional
    %variables, e.g. stepsize, timestep, local work group size, etc...

    global I_Mesh I_TI I_BalanceLaws I_Tech I_RunOps
    
    if strcmp(I_RunOps('periodic'), 'USE_PERIODIC')
        % the nodes at the right boundary are not included
        DX = double(I_Mesh('XMAX') - I_Mesh('XMIN')) / double(I_Mesh('NODES_X'));
        DY = double(I_Mesh('YMAX') - I_Mesh('YMIN')) / double(I_Mesh('NODES_Y'));
        DZ = double(I_Mesh('ZMAX') - I_Mesh('ZMIN')) / double(I_Mesh('NODES_Z'));
    else
        DX = double(I_Mesh('XMAX') - I_Mesh('XMIN')) / (double(I_Mesh('NODES_X'))-1);
        DY = double(I_Mesh('YMAX') - I_Mesh('YMIN')) / (double(I_Mesh('NODES_Y'))-1);
        DZ = double(I_Mesh('ZMAX') - I_Mesh('ZMIN')) / (double(I_Mesh('NODES_Z'))-1);
    end

    I_Mesh('DX') = DX; I_Mesh('DY') = DY; I_Mesh('DZ') = DZ;

    I_BalanceLaws('g_range') = uint32([I_Mesh('NODES_X'), I_Mesh('NODES_Y'), I_Mesh('NODES_Z')]);
    I_BalanceLaws('l_range') = uint32([0]);

    num_nodes = I_Mesh('NODES_X')*I_Mesh('NODES_Y')*I_Mesh('NODES_Z');
    
    if strcmp(I_Tech('REAL'),'float')
        field_u1 = single(zeros(1, num_nodes*I_BalanceLaws('NUM_TOTAL_VARS')));
        field_u2 = single(zeros(1, num_nodes*I_BalanceLaws('NUM_TOTAL_VARS')));
    else
        field_u1 = double(zeros(1, num_nodes*I_BalanceLaws('NUM_TOTAL_VARS')));
        field_u2 = double(zeros(1, num_nodes*I_BalanceLaws('NUM_TOTAL_VARS')));
    end
        

    x = linspace(I_Mesh('XMIN'), I_Mesh('XMAX'), I_Mesh('NODES_X'));
    y = linspace(I_Mesh('YMIN'), I_Mesh('YMAX'), I_Mesh('NODES_Y'));
    z = linspace(I_Mesh('ZMIN'), I_Mesh('ZMAX'), I_Mesh('NODES_Z'));
    
    kernel_path_list = {};
    if I_RunOps('order') == 2
        if strcmp(I_RunOps('operator_form'), 'extended')
            kernel_path_list = [kernel_path_list, {'../include/2ndOrderExtended.h'}];
        else % classical operators
            kernel_path_list = [kernel_path_list, {'../include/2ndOrder.h'}];
        end
    elseif I_RunOps('order') == 4
        if strcmp(I_RunOps('operator_form'), 'extended')
            kernel_path_list = [kernel_path_list, {'../include/4thOrderExtended.h'}];
        else % classical operators
            kernel_path_list = [kernel_path_list, {'../include/4thOrder.h'}];
        end
    elseif I_RunOps('order') == 6
        if strcmp(I_RunOps('operator_form'), 'extended')
            kernel_path_list = [kernel_path_list, {'../include/6thOrderExtended.h'}];
        else % classical operators
            kernel_path_list = [kernel_path_list, {'../include/6thOrder.h'}];
        end
    else
        fprintf('Specify order \n')
    end
    kernel_path_list = [kernel_path_list, {sprintf('../include_defines/%s_defines.h', I_RunOps('conservation_laws'))}];
    kernel_path_list = [kernel_path_list, {'../include/utils.h'}];
    kernel_path_list = [kernel_path_list, {'../kernel/SBP_operator.cl'}];
    kernel_path_list = [kernel_path_list, {sprintf('../include_testcases/%s_%s.h', I_RunOps('conservation_laws'), I_RunOps('testcase'))}];
    kernel_path_list = [kernel_path_list, {sprintf('../include_physics/%s.h', I_RunOps('conservation_laws'))}];
    kernel_path_list = [kernel_path_list, {'../kernel/kernel_init.cl'}];
    
    settings_tech = generate_settings(I_Tech, {'REAL'; 'REAL4'; 'optimizations'});
    settings_mesh = generate_settings(I_Mesh, {'DX'; 'DY'; 'DZ';...
                                               'NODES_X'; 'NODES_Y'; 'NODES_Z';...
                                               'XMIN'; 'XMAX'; 'YMIN'; 'YMAX'; 'ZMAX'; 'ZMIN'});
    settings_runops = generate_settings(I_RunOps, {'periodic'});

    settings = strcat(settings_tech, settings_mesh, settings_runops);
    
    [~] = cl_run_kernel(I_Tech('device'), kernel_path_list, settings);
    
    cl_run_kernel(I_Tech('device'), 'init', I_BalanceLaws('g_range'), I_BalanceLaws('l_range'), field_u1, 0);

    kernel_path_list = {};
    % Include header file containing the coefficients of the respective order
    if I_RunOps('order') == 2
        if strcmp(I_RunOps('operator_form'), 'extended')
            kernel_path_list = [kernel_path_list, {'../include/2ndOrderExtended.h'}];
        else % classical operators
            kernel_path_list = [kernel_path_list, {'../include/2ndOrder.h'}];
        end
    elseif I_RunOps('order') == 4
        if strcmp(I_RunOps('operator_form'), 'extended')
            kernel_path_list = [kernel_path_list, {'../include/4thOrderExtended.h'}];
        else % classical operators
            kernel_path_list = [kernel_path_list, {'../include/4thOrder.h'}];
        end
    elseif I_RunOps('order') == 6
        if strcmp(I_RunOps('operator_form'), 'extended')
            kernel_path_list = [kernel_path_list, {'../include/6thOrderExtended.h'}];
        else % classical operators
            kernel_path_list = [kernel_path_list, {'../include/6thOrder.h'}];
        end
    else
        fprintf('Specify order \n')
    end
    kernel_path_list = [kernel_path_list, {sprintf('../include_defines/%s_defines.h', I_RunOps('conservation_laws'))}];
    kernel_path_list = [kernel_path_list, {'../include/utils.h'}];
    kernel_path_list = [kernel_path_list, {'../kernel/SBP_operator.cl'}];
    kernel_path_list = [kernel_path_list, {sprintf('../include_testcases/%s_%s.h', I_RunOps('conservation_laws'), I_RunOps('testcase'))}];
    kernel_path_list = [kernel_path_list, {sprintf('../include_physics/%s.h', I_RunOps('conservation_laws'))}];
    kernel_path_list = [kernel_path_list, {'../kernel/kernel_time_integrator.cl'}];

    settings_tech = generate_settings(I_Tech, {'REAL'; 'REAL4'; 'optimizations'});
    settings_mesh = generate_settings(I_Mesh, {'DX'; 'DY'; 'DZ';...
                                               'NODES_X'; 'NODES_Y'; 'NODES_Z';...
                                               'XMIN'; 'XMAX'; 'YMIN'; 'YMAX'; 'ZMAX'; 'ZMIN'});
    settings_time_integration = generate_settings(I_TI, {'DT'});
    settings_runops = generate_settings(I_RunOps, {'periodic'});

    settings = strcat(settings_tech, settings_mesh, settings_time_integration, settings_runops);

    % Compile kernel
    [~] = cl_run_kernel(I_Tech('device'), kernel_path_list, settings);

end
