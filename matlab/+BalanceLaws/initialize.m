%This project is licensed under the terms of the Creative Commons CC BY-NC-ND 4.0 license.

function [field_u1, field_u2] = initialize()

    %Initialize the magnetic field, the velocity field and the density field
    %according to the specified testcase. Also calculates and sets additional
    %variables, e.g. stepsize, timestep, local work group size, etc...

    global I_Mesh I_TI I_BalanceLaws I_Tech I_RunOps

    if strcmp(I_RunOps('periodic_x'), 'USE_PERIODIC_X')
        % the nodes at the right boundary are not included
        DX = double(I_Mesh('XMAX') - I_Mesh('XMIN')) / double(I_Mesh('NODES_X'));
    else
        DX = double(I_Mesh('XMAX') - I_Mesh('XMIN')) / (double(I_Mesh('NODES_X'))-1);
    end
    if strcmp(I_RunOps('periodic_y'), 'USE_PERIODIC_Y')
        % the nodes at the right boundary are not included
        DY = double(I_Mesh('YMAX') - I_Mesh('YMIN')) / double(I_Mesh('NODES_Y'));
    else
        DY = double(I_Mesh('YMAX') - I_Mesh('YMIN')) / (double(I_Mesh('NODES_Y'))-1);
    end
    if strcmp(I_RunOps('periodic_z'), 'USE_PERIODIC_Z')
        % the nodes at the right boundary are not included
        DZ = double(I_Mesh('ZMAX') - I_Mesh('ZMIN')) / double(I_Mesh('NODES_Z'));
    else
        DZ = double(I_Mesh('ZMAX') - I_Mesh('ZMIN')) / (double(I_Mesh('NODES_Z'))-1);
    end

    I_Mesh('DX') = DX; I_Mesh('DY') = DY; I_Mesh('DZ') = DZ;

    I_BalanceLaws('g_range') = uint32([I_Mesh('NODES_X'), I_Mesh('NODES_Y'), I_Mesh('NODES_Z')]);
    I_BalanceLaws('l_range') = uint32([0]);

    num_nodes = I_Mesh('NODES_X')*I_Mesh('NODES_Y')*I_Mesh('NODES_Z');
%%
    % Depending on the device type choose the local work group size. If the
    % device is a GPU the highest availble work group size is selected.
    % For CPUs the optimal work group size was determined
    % through testing. Especially for high-end Intel CPUs a higher value can decrease runtime.
    [~, dev_type, ~, ~, lw_size, ~] = cl_get_devices;
    type = dev_type(I_Tech('device'));
    if (strcmp(type{1},'CPU'))
        group_size = 16;
    else
        if (num_nodes > lw_size(I_Tech('device')))
            group_size = lw_size(I_Tech('device'));
        else
            group_size = 2^floor(log(double(num_nodes)) / log(2));
        end
    end
    group_size = double(group_size);
    num_nodes_pad = ceil(double(num_nodes)/group_size)*group_size;
    num_groups = ceil(num_nodes_pad/group_size);

    I_Tech('NUM_NODES_PAD') = uint32(num_nodes_pad);
    I_Tech('num_groups') = num_groups;
    I_Tech('W_SIZE') = uint16(group_size);

    % Global and local range for dot product and norm
    I_Tech('g_range') = uint32([I_Tech('NUM_NODES_PAD'), 1, 1]);
    I_Tech('l_range') = uint32([I_Tech('W_SIZE'), 1, 1]);
%%
    if strcmp(I_Tech('REAL'),'float')
        field_u1 = single(zeros(1, I_Tech('NUM_NODES_PAD')*I_BalanceLaws('NUM_TOTAL_VARS')));
        field_u2 = single(zeros(1, I_Tech('NUM_NODES_PAD')*I_BalanceLaws('NUM_TOTAL_VARS')));
    else
        field_u1 = double(zeros(1, I_Tech('NUM_NODES_PAD')*I_BalanceLaws('NUM_TOTAL_VARS')));
        field_u2 = double(zeros(1, I_Tech('NUM_NODES_PAD')*I_BalanceLaws('NUM_TOTAL_VARS')));
    end

%%
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

    settings_tech = generate_settings(I_Tech, {'REAL'; 'REAL4'; 'optimizations'; 'memory_layout'; 'NUM_NODES_PAD'});
    settings_mesh = generate_settings(I_Mesh, {'DX'; 'DY'; 'DZ';...
                                               'NODES_X'; 'NODES_Y'; 'NODES_Z';...
                                               'XMIN'; 'XMAX'; 'YMIN'; 'YMAX'; 'ZMAX'; 'ZMIN'});
    settings_runops = generate_settings(I_RunOps, {'periodic_x'; 'periodic_y'; 'periodic_z'});
    settings_time_integration = generate_settings(I_TI, {'DT'});

    settings = strcat(settings_tech, settings_mesh, settings_runops, settings_time_integration);

    [~] = cl_run_kernel(I_Tech('device'), kernel_path_list, settings);

    cl_run_kernel(I_Tech('device'), 'init', I_BalanceLaws('g_range'), I_BalanceLaws('l_range'), field_u1, 0);
%%
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
    kernel_path_list = [kernel_path_list, {'../kernel/kernel_norm.cl'}];

    settings_tech = generate_settings(I_Tech, {'REAL'; 'REAL4'; 'W_SIZE'; 'optimizations'; 'memory_layout'; 'NUM_NODES_PAD'});
    settings_mesh = generate_settings(I_Mesh, {'DX'; 'DY'; 'DZ';...
                                               'NODES_X'; 'NODES_Y'; 'NODES_Z';...
                                               'XMIN'; 'XMAX'; 'YMIN'; 'YMAX'; 'ZMAX'; 'ZMIN'});
    settings_time_integration = generate_settings(I_TI, {'DT'});
    settings_runops = generate_settings(I_RunOps, {'periodic_x'; 'periodic_y'; 'periodic_z'});

    settings = strcat(settings_tech, settings_mesh, settings_time_integration, settings_runops);

    % Compile kernel
    [~] = cl_run_kernel(I_Tech('device'), kernel_path_list, settings);

end
