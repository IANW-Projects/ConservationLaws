%This project is licensed under the terms of the Creative Commons CC BY-NC-ND 4.0 license.

function [] = compute_numerical_solution(field_u1, field_u2)

% Takes the initialized fields and advances the solution in time

global I_Mesh I_TI I_BalanceLaws I_Tech I_RunOps I_Results

if strcmp(I_Tech('REAL'),'float')
    current_time = single(zeros(2));
    components = single(zeros(2));
    norm_output = single(zeros(1, I_Tech('num_groups')));
    Lerror = single(zeros(I_BalanceLaws('NUM_CONSERVED_VARS'), I_TI('num_steps')+1));
else
    current_time = double(zeros(2));
    components = double(zeros(2));
    norm_output = double(zeros(1, I_Tech('num_groups')));
    Lerror = double(zeros(I_BalanceLaws('NUM_CONSERVED_VARS'), I_TI('num_steps')+1));
end

switch I_TI('time_integrator')

    case 'CarpenterKennedy2N54'
        RK_Step_a = {'CarpenterKennedy2N54_1a', 'CarpenterKennedy2N54_2a', ...
                     'CarpenterKennedy2N54_3a', 'CarpenterKennedy2N54_4a', ...
                     'CarpenterKennedy2N54_5a'};
        RK_Step_b = {'CarpenterKennedy2N54_1b', 'CarpenterKennedy2N54_2b', ...
                     'CarpenterKennedy2N54_3b', 'CarpenterKennedy2N54_4b', ...
                     'CarpenterKennedy2N54_5b'};
        RK_aux_vars = repmat({'calc_auxiliary_vars_1_2_args'}, size(RK_Step_a));

        tmp = [RK_aux_vars; RK_Step_a; RK_Step_b];
        RK_Step = tmp(:)';
        time_integrator_num_fields = 2;

    case 'ToulorgeDesmet2N84F'

        RK_Step_a = {'ToulorgeDesmet2N84F_1a', 'ToulorgeDesmet2N84F_2a', ...
                     'ToulorgeDesmet2N84F_3a', 'ToulorgeDesmet2N84F_4a', ...
                     'ToulorgeDesmet2N84F_5a', 'ToulorgeDesmet2N84F_6a', ...
                     'ToulorgeDesmet2N84F_7a', 'ToulorgeDesmet2N84F_8a'};
        RK_Step_b = {'ToulorgeDesmet2N84F_1b', 'ToulorgeDesmet2N84F_2b', ...
                     'ToulorgeDesmet2N84F_3b', 'ToulorgeDesmet2N84F_4b', ...
                     'ToulorgeDesmet2N84F_5b', 'ToulorgeDesmet2N84F_6b', ...
                     'ToulorgeDesmet2N84F_7b', 'ToulorgeDesmet2N84F_8b'};
        RK_aux_vars = repmat({'calc_auxiliary_vars_1_2_args'}, size(RK_Step_a));
        tmp = [RK_aux_vars; RK_Step_a; RK_Step_b];
        RK_Step = tmp(:)';
        time_integrator_num_fields = 2;

    case 'KennedyCarpenterLewis2R54C'
        RK_Step_a = {'KennedyCarpenterLewis2R54C_1a', 'KennedyCarpenterLewis2R54C_2a', ...
                     'KennedyCarpenterLewis2R54C_3a', 'KennedyCarpenterLewis2R54C_4a', ...
                     'KennedyCarpenterLewis2R54C_5'};
        RK_Step_b = {'KennedyCarpenterLewis2R54C_1b', 'KennedyCarpenterLewis2R54C_2b', ...
                     'KennedyCarpenterLewis2R54C_3b', 'KennedyCarpenterLewis2R54C_4b', ...
                     'calc_time_3_args'};

        RK_aux_vars = {'calc_auxiliary_vars_1_3_args', 'calc_auxiliary_vars_1_3_args', 'calc_auxiliary_vars_2_3_args',...
                       'calc_auxiliary_vars_1_3_args', 'calc_auxiliary_vars_2_3_args'};

        tmp = [RK_aux_vars; RK_Step_a; RK_Step_b];
        RK_Step = tmp(:)';
        time_integrator_num_fields = 3;

    case 'CalvoFrancoRandez2R64'

        RK_Step_a = {'CalvoFrancoRandez2R64_1a', 'CalvoFrancoRandez2R64_2a', ...
                     'CalvoFrancoRandez2R64_3a', 'CalvoFrancoRandez2R64_4a', ...
                     'CalvoFrancoRandez2R64_5a', 'CalvoFrancoRandez2R64_6'};
        RK_Step_b = {'CalvoFrancoRandez2R64_1b', 'CalvoFrancoRandez2R64_2b', ...
                     'CalvoFrancoRandez2R64_3b', 'CalvoFrancoRandez2R64_4b', ...
                     'CalvoFrancoRandez2R64_5b', 'calc_time_3_args'};

        RK_aux_vars = {'calc_auxiliary_vars_1_3_args', 'calc_auxiliary_vars_2_3_args', 'calc_auxiliary_vars_2_3_args',...
                       'calc_auxiliary_vars_2_3_args', 'calc_auxiliary_vars_2_3_args', 'calc_auxiliary_vars_2_3_args'};
        tmp = [RK_aux_vars; RK_Step_a; RK_Step_b];
        RK_Step = tmp(:)';
        time_integrator_num_fields = 3;

    case 'SSPRK33'
        RK_Step = [];
        RK_Step{end + 1} = 'calc_auxiliary_vars_1_3_args';
        RK_Step{end + 1} = 'SSPRK33_1';
        RK_Step{end + 1} = 'SSPRK33_2a';
        RK_Step{end + 1} = 'calc_auxiliary_vars_2_3_args';
        RK_Step{end + 1} = 'SSPRK33_2b';
        RK_Step{end + 1} = 'SSPRK33_3a';
        RK_Step{end + 1} = 'calc_auxiliary_vars_3_3_args';
        RK_Step{end + 1} = 'SSPRK33_3b';
        RK_Step{end + 1} = 'calc_time_3_args';
        time_integrator_num_fields = 3;

    case 'SSPRK104'
        RK_Step = [];
        RK_Step{end + 1} = 'calc_auxiliary_vars_1_3_args';
        RK_Step{end + 1} = 'SSPRK104_01';
        RK_Step{end + 1} = 'SSPRK104_02a';
        RK_Step{end + 1} = 'calc_auxiliary_vars_2_3_args';
        RK_Step{end + 1} = 'SSPRK104_02b';
        RK_Step{end + 1} = 'SSPRK104_03a';
        RK_Step{end + 1} = 'calc_auxiliary_vars_3_3_args';
        RK_Step{end + 1} = 'SSPRK104_03b';
        RK_Step{end + 1} = 'SSPRK104_04a';
        RK_Step{end + 1} = 'calc_auxiliary_vars_2_3_args';
        RK_Step{end + 1} = 'SSPRK104_04b';
        RK_Step{end + 1} = 'SSPRK104_05a';
        RK_Step{end + 1} = 'calc_auxiliary_vars_3_3_args';
        RK_Step{end + 1} = 'SSPRK104_05b';
        RK_Step{end + 1} = 'SSPRK104_06';
        RK_Step{end + 1} = 'calc_auxiliary_vars_2_3_args';
        RK_Step{end + 1} = 'SSPRK104_07';
        RK_Step{end + 1} = 'SSPRK104_08a';
        RK_Step{end + 1} = 'calc_auxiliary_vars_1_3_args';
        RK_Step{end + 1} = 'SSPRK104_08b';
        RK_Step{end + 1} = 'SSPRK104_09a';
        RK_Step{end + 1} = 'calc_auxiliary_vars_2_3_args';
        RK_Step{end + 1} = 'SSPRK104_09b';
        RK_Step{end + 1} = 'SSPRK104_10a';
        RK_Step{end + 1} = 'calc_auxiliary_vars_1_3_args';
        RK_Step{end + 1} = 'SSPRK104_10b';
        RK_Step{end + 1} = 'SSPRK104_11a';
        RK_Step{end + 1} = 'calc_auxiliary_vars_2_3_args';
        RK_Step{end + 1} = 'SSPRK104_11b';
        time_integrator_num_fields = 3;
end

num_steps_run = I_TI('num_steps');
kernel_runtime = 0;
RK_block_size = 15000;

%%
tic
switch time_integrator_num_fields
    case 2
        if I_RunOps('save_integrals_over_time')
            for step = 1:I_TI('num_steps')
                for comp=0:I_BalanceLaws('NUM_CONSERVED_VARS')-1
                    components(1) = comp;
                    cl_run_kernel(I_Tech('device'), 'analytical_u', I_BalanceLaws('g_range'), I_BalanceLaws('l_range'), field_u2, current_time, 0);
                    norm_output(:) = 0;
                    if strcmp(I_RunOps('norm'),'L2')
                        cl_run_kernel(I_Tech('device'), 'norm2_diff', I_Tech('g_range'), I_Tech('l_range'), field_u1, field_u2, norm_output, components, 0);
                        Lerror(comp + 1, step) = sqrt(sum(norm_output));
                    elseif strcmp(I_RunOps('norm'),'LInf')
                        cl_run_kernel(I_Tech('device'), 'norm_infty_diff', I_Tech('g_range'), I_Tech('l_range'), field_u1, field_u2, norm_output, components, 0);
                        Lerror(comp + 1, step) = max(norm_output);
                    end
                end

                t = cl_run_kernel(I_Tech('device'), RK_Step, I_BalanceLaws('g_range'), I_BalanceLaws('l_range'), ...
                                      field_u1, field_u2, ...
                                      current_time, 0);
            end

        else
            while num_steps_run > RK_block_size
                kernel_list = repmat(RK_Step, 1, RK_block_size);

                t = cl_run_kernel(I_Tech('device'), kernel_list, I_BalanceLaws('g_range'), I_BalanceLaws('l_range'), ...
                                  field_u1, field_u2, ...
                                  current_time, 0);
                kernel_runtime =  kernel_runtime + t;
                num_steps_run = num_steps_run - RK_block_size;
            end
            if num_steps_run > 0
                kernel_list = repmat(RK_Step, 1, num_steps_run);
                t = cl_run_kernel(I_Tech('device'), kernel_list, I_BalanceLaws('g_range'), I_BalanceLaws('l_range'), ...
                                              field_u1, field_u2, ...
                                              current_time, 0);
                kernel_runtime = kernel_runtime + t;
            end
        end
    case 3
        num_nodes = I_Mesh('NODES_X')*I_Mesh('NODES_Y')*I_Mesh('NODES_Z');
        if strcmp(I_Tech('REAL'),'float')
            field_u3 = single(zeros(1, I_Tech('num_nodes_pad')*I_BalanceLaws('NUM_TOTAL_VARS')));
        else
            field_u3 = double(zeros(1, I_Tech('num_nodes_pad')*I_BalanceLaws('NUM_TOTAL_VARS')));
        end

        if I_RunOps('save_integrals_over_time')
            for step = 1:I_TI('num_steps')
                for comp=0:I_BalanceLaws('NUM_CONSERVED_VARS')-1
                    components(1) = comp;
                    cl_run_kernel(I_Tech('device'), 'analytical_u', I_BalanceLaws('g_range'), I_BalanceLaws('l_range'), field_u2, current_time, 0);
                    norm_output(:) = 0;
                    if strcmp(I_RunOps('norm'),'L2')
                        cl_run_kernel(I_Tech('device'), 'norm2_diff', I_Tech('g_range'), I_Tech('l_range'), field_u1, field_u2, norm_output, components, 0);
                        Lerror(comp + 1, step) = sqrt(sum(norm_output));
                    elseif strcmp(I_RunOps('norm'),'LInf')
                        cl_run_kernel(I_Tech('device'), 'norm_infty_diff', I_Tech('g_range'), I_Tech('l_range'), field_u1, field_u2, norm_output, components, 0);
                        Lerror(comp + 1, step) = max(norm_output);
                    end
                end

                t = cl_run_kernel(I_Tech('device'), RK_Step, I_BalanceLaws('g_range'), I_BalanceLaws('l_range'), ...
                                      field_u1, field_u2, field_u3, ...
                                      current_time, 0);
            end

        else
            while num_steps_run > RK_block_size
                kernel_list = repmat(RK_Step, 1, RK_block_size);

                t = cl_run_kernel(I_Tech('device'), kernel_list, I_BalanceLaws('g_range'), I_BalanceLaws('l_range'), ...
                                  field_u1, field_u2, field_u3, ...
                                  current_time, 0);
                kernel_runtime =  kernel_runtime + t;
                num_steps_run = num_steps_run - RK_block_size;
            end
            if num_steps_run > 0
                kernel_list = repmat(RK_Step, 1, num_steps_run);
                t = cl_run_kernel(I_Tech('device'), kernel_list, I_BalanceLaws('g_range'), I_BalanceLaws('l_range'), ...
                                 field_u1, field_u2, field_u3, ...
                                 current_time, 0);
                kernel_runtime = kernel_runtime + t;
            end
        end
end
runtime = toc;

% save runtime and kernel_runtime
I_Results('runtime') = runtime;
I_Results('kernel_runtime') = kernel_runtime;

% save fields
if I_RunOps('save_fields')
    I_Results('field_u') = field_u1;
end

%Calculate analytical solution and Lerror
current_time(1) = I_TI('final_time');
cl_run_kernel(I_Tech('device'), 'analytical_u', I_BalanceLaws('g_range'), I_BalanceLaws('l_range'), field_u2, current_time, 0);

abs_err = zeros(I_BalanceLaws('NUM_CONSERVED_VARS'),1);
rel_err = zeros(I_BalanceLaws('NUM_CONSERVED_VARS'),1);

for comp=0:I_BalanceLaws('NUM_CONSERVED_VARS')-1
	components(1) = comp;
    norm_output(:) = 0;
    if strcmp(I_RunOps('norm'),'L2')
        cl_run_kernel(I_Tech('device'), 'norm2_diff', I_Tech('g_range'), I_Tech('l_range'), field_u1, field_u2, norm_output, components, 0);
        abs_err(comp + 1) = sqrt(sum(norm_output));

        norm_output(:) = 0;
        cl_run_kernel(I_Tech('device'), 'norm2', I_Tech('g_range'), I_Tech('l_range'), field_u2, norm_output, components, 0);
        rel_err(comp + 1) = abs_err(comp + 1) / sqrt(sum(norm_output));
    elseif strcmp(I_RunOps('norm'),'LInf')
        cl_run_kernel(I_Tech('device'), 'norm_infty_diff', I_Tech('g_range'), I_Tech('l_range'), field_u1, field_u2, norm_output, components, 0);
        abs_err(comp + 1) = max(norm_output);

        norm_output(:) = 0;
        cl_run_kernel(I_Tech('device'), 'norm_infty', I_Tech('g_range'), I_Tech('l_range'), field_u2, norm_output, components, 0);
        rel_err(comp + 1) = abs_err(comp + 1) / max(norm_output);
    end
end

I_Results('abs_err') = abs_err;
I_Results('rel_err') = rel_err;

if I_RunOps('save_integrals_over_time')
    Lerror(:,I_TI('num_steps')+1) = I_Results('abs_err');
    I_Results('Lerror_over_time') = Lerror;
    I_Results('time') = linspace(0, I_TI('final_time'), size(Lerror,2));
end

end
