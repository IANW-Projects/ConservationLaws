import numpy as np
import sys
import math

prepare_vars()


N_fac = 1

I_Mesh['NODES_X'] = int(N_fac * 100)
I_Mesh['NODES_Y'] = int(N_fac * 100)
I_Mesh['NODES_Z'] = int(N_fac * 100)

I_Mesh['XMIN'] = 0.0
I_Mesh['YMIN'] = 0.0
I_Mesh['ZMIN'] = 0.0

I_Mesh['XMAX'] = 1.0
I_Mesh['YMAX'] = 1.0
I_Mesh['ZMAX'] = 1.0

I_TI['final_time'] = 0.1
I_TI['cfl'] = 0.2
k = 10

dt = I_TI['cfl'] * 2.0 / float(I_Mesh['NODES_Y']) # this has to be estimated
num_steps = math.ceil(I_TI['final_time']/dt)
dt = I_TI['final_time'] / num_steps

I_TI['time_integrator'] = 'CarpenterKennedy2N54'

I_TI['DT'] = dt
I_TI['num_steps'] = num_steps

# I_Tech['device'] = 1

I_Tech['REAL'] = 'float' #'double'
I_Tech['REAL4'] = I_Tech['REAL'] + "4"
I_Tech['memory_layout'] = 'USE_ARRAY_OF_STRUCTURES' #'USE_STRUCTURE_OF_ARRAYS'


I_BalanceLaws['NUM_CONSERVED_VARS'] = 8
I_BalanceLaws['NUM_AUXILIARY_VARS'] = 4
I_BalanceLaws['NUM_TOTAL_VARS'] = I_BalanceLaws['NUM_CONSERVED_VARS'] + I_BalanceLaws['NUM_AUXILIARY_VARS']

# Compiler based optimizations

if I_Tech['REAL'] == 'float':
    I_Tech['optimizations'] = ' -cl-mad-enable -cl-no-signed-zeros -cl-finite-math-only -cl-single-precision-constant'
else:
    I_Tech['optimizations'] = ' -cl-mad-enable -cl-no-signed-zeros -cl-finite-math-only'

I_RunOps['periodicx'] = 'none'# 'USE_PERIODIC'
I_RunOps['periodicy'] = 'none'
I_RunOps['periodicz'] = 'none'

I_RunOps['vr'] = 'none'#'USE_VR_KUSANO'

I_RunOps['order'] = 2
I_RunOps['conservation_laws'] = 'ideal_MHD'
I_RunOps['testcase'] = 'far_dipole'
I_RunOps['plot_numerical_solution'] = 'z'
I_RunOps['save_integrals_over_time'] = False
I_RunOps['norm'] = 'LInf'

## Initialize

field_u1, field_u2 = initialize()

if len(sys.argv) > 1:
    field_u1 = reload_all_variables(field_u1, sys.argv[1], ['rho', 'px', 'py', 'pz', 'E', 'Bx', 'By', 'Bz'])
    print(f"reloaded the conserved variables out of {sys.argv[1]}")

print('Testcase: ' + I_RunOps['testcase'] + '\nOrder: ' + str(I_RunOps['order']) + '\nTime integrator: ' + \
        I_TI['time_integrator'] +' \nDT: ' + str(I_TI['DT']) + '\nN_STEPS: ' +str(I_TI['num_steps']) + " FINAL_TIME: " + \
        str(I_TI['final_time']) + '\nDX: ' + str(I_Mesh['DX']) + ' NODES_X: ' +str(I_Mesh['NODES_X'])+'\nDY: '+ str(I_Mesh['DY']) + \
        ' NODES_Y: ' + str(I_Mesh['NODES_Y']) + '\nDZ: ' + str(I_Mesh['DZ'])+  ' NODES_Z: ' + str(I_Mesh['NODES_Z']) + '\nREAL: ' + \
        str(I_Tech['REAL']))

print(str(k) + " snapshots")
ct = False

for i in range(k):
    ct = compute_numerical_solution(field_u1, field_u2, ct)
    save_all_variables(field_u1, "results/output" + str(i) + ".hdf5", ['rho', 'px', 'py', 'pz', 'E', 'Bx', 'By', 'Bz'])
    print(i)

#rel_err = I_Results['rel_err']
#for comp in range(I_BalanceLaws['NUM_CONSERVED_VARS']):
 #   print("Relative Error of Field Component " + str(comp) +": " + str(100*rel_err[comp]))


#field_u1_reshaped = np.reshape(field_u1, (I_Tech['NUM_NODES_PAD'], I_BalanceLaws['NUM_TOTAL_VARS']))



#plot_fieldlines_B(field_u1_reshaped, I_RunOps['plot_numerical_solution'], I_Mesh['NODES_X'], I_Mesh['NODES_Y'], I_Mesh['NODES_Z'], 'Numerical Solution', Field_Bx, Field_By)
#plot_2D(field_u1_reshaped, I_RunOps['plot_numerical_solution'], I_Mesh['NODES_X'], I_Mesh['NODES_Y'], I_Mesh['NODES_Z'], 'Numerical Solution', Field_By,Field_By)