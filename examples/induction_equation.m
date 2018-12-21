clc
clear
close all

addpath('../matlab')

BalanceLaws.prepare_vars();
global I_Mesh I_TI I_BalanceLaws I_Tech I_RunOps I_Results

N = uint32(40);
I_Mesh('NODES_X') = N; I_Mesh('NODES_Y') = N; I_Mesh('NODES_Z') = N;

I_Mesh('XMIN') = 0; I_Mesh('XMAX') = 4*pi/3;
I_Mesh('YMIN') = 0; I_Mesh('YMAX') = 4*pi/3;
I_Mesh('ZMIN') = 0; I_Mesh('ZMAX') = 4*pi/3;

I_TI('final_time') = 1;
I_TI('cfl') = 0.95/double(N); %nonperiodic boundaries: 0.95/double(N)
                              %periodic boundaries: 2/double(N)

%Set dt
dt = 1.0405827263267431e-03; %periodic:2.1367521367521370e-03 ; nonperiodic: 1.0405827263267431e-03
num_steps = ceil(I_TI('final_time')/dt);
dt = I_TI('final_time') / num_steps;

I_TI('DT') = dt;
I_TI('num_steps') = num_steps;
%Time integrator: SSPRK33 SSPRK104 CalvoFrancoRandez2R64
%KennedyCarpenterLewis2R54C CarpenterKennedy2N54 ToulorgeDesmet2N84F
I_TI('time_integrator') = 'CarpenterKennedy2N54';

I_Tech('device') = 1;
I_Tech('REAL') = 'double'; % float
I_Tech('REAL4') = sprintf('%s4',I_Tech('REAL')); %Vector datatype

% Use as kernel defines to keep consistency with header files?
I_BalanceLaws('NUM_CONSERVED_VARS') = 3;
I_BalanceLaws('NUM_AUXILIARY_VARS') = 7;
I_BalanceLaws('NUM_TOTAL_VARS') = I_BalanceLaws('NUM_CONSERVED_VARS') + I_BalanceLaws('NUM_AUXILIARY_VARS');

%Compiler based optimizations
if strcmp(I_Tech('REAL'),'float')
    I_Tech('optimizations') = ' -cl-mad-enable -cl-no-signed-zeros -cl-finite-math-only -cl-single-precision-constant';
else
    I_Tech('optimizations') = ' -cl-mad-enable -cl-no-signed-zeros -cl-finite-math-only';
end

I_RunOps('periodic') = 'NONE'; % 'NONE', 'USE_PERIODIC'; must be set to 'USE_PERIODIC'
                                       % if periodic boundary conditions should be used

I_RunOps('order') = 4; I_RunOps('operator_form') = 'classical'; % order: 2, 4, 6; operator_form: classical, extended
I_RunOps('conservation_laws') = 'induction_equation';
I_RunOps('testcase') = 'hall_periodic';
I_RunOps('plot_numerical_solution') = 'z';
I_RunOps('save_integrals_over_time') = false;
% Choose between L2 and LInfinity norm for error calculation
I_RunOps('norm') = 'LInf'; % L2, LInf
%% Initialize variables
[field_u1, field_u2] = BalanceLaws.initialize();

fprintf('Testcase: %s \nOrder: %d \nTime integrator: %s\nDT: %.16e   N_STEPS: %5d   FINAL_TIME: %.16e\nDX: %.16e   NODES_X: %5d\nDY: %.16e   NODES_Y: %5d\nDZ: %.16e   NODES_Z: %5d \nREAL: %s\n\n',...
        I_RunOps('testcase'), I_RunOps('order'), I_TI('time_integrator'), I_TI('DT'), I_TI('num_steps'), I_TI('final_time'), I_Mesh('DX'), I_Mesh('NODES_X'), I_Mesh('DY'), I_Mesh('NODES_Y'), I_Mesh('DZ'), I_Mesh('NODES_Z'), I_Tech('REAL'));
%% Compute numerical solution
BalanceLaws.compute_numerical_solution(field_u1, field_u2);

if strcmp(I_RunOps('periodic'), 'USE_PERIODIC')
    rel_err = I_Results('rel_err');
    for comp=0:I_BalanceLaws('NUM_CONSERVED_VARS') - 1
        fprintf('Relative Error of Field Component %d: %.15f %%\n', comp, 100*rel_err(comp + 1))
    end
end

%% Plot numerical solution
field_u1_reshaped = reshape(field_u1, [I_BalanceLaws('NUM_TOTAL_VARS'), I_Tech('num_nodes_pad')]);

%Optional plots
if ismember(lower(char(I_RunOps('plot_numerical_solution'))),{'x','y','z','xy', 'xz', 'yz', 'xyz'})
    plot_2D(field_u1_reshaped, I_RunOps('plot_numerical_solution'),...
        I_Mesh('NODES_X'), I_Mesh('NODES_Y'), I_Mesh('NODES_Z'), 'Numerical Solution', 1, 3);
end
