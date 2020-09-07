import numpy as np
import sys
import math



prepare_vars()


N_fac = 1

I_Mesh['NODES_X'] = int(N_fac * 58)
I_Mesh['NODES_Y'] = int(N_fac * 100)
I_Mesh['NODES_Z'] = int(N_fac * 12)

I_Mesh['XMIN'] = 0.0
I_Mesh['YMIN'] = 0.0
I_Mesh['ZMIN'] = 0.0

I_Mesh['XMAX'] = 2.0/math.sqrt(3)
I_Mesh['YMAX'] = 2.0
I_Mesh['ZMAX'] = 1.0

I_TI['final_time'] = 5
I_TI['cfl'] = 0.3

dt = I_TI['cfl'] * 2.0 / float(I_Mesh['NODES_Y'])
num_steps = math.ceil(I_TI['final_time']/dt)
dt = I_TI['final_time'] / num_steps

I_TI['time_integrator'] = 'CarpenterKennedy2N54'

I_TI['DT'] = dt
I_TI['num_steps'] = num_steps

# I_Tech['device'] = 1

I_Tech['REAL'] = 'double'
I_Tech['REAL4'] = I_Tech['REAL'] + "4"
I_Tech['memory_layout'] = 'USE_ARRAY_OF_STRUCTURES'#'USE_STRUCTURE_OF_ARRAYS'


I_BalanceLaws['NUM_CONSERVED_VARS'] = 8
I_BalanceLaws['NUM_AUXILIARY_VARS'] = 4
I_BalanceLaws['NUM_TOTAL_VARS'] = I_BalanceLaws['NUM_CONSERVED_VARS'] + I_BalanceLaws['NUM_AUXILIARY_VARS']

# Compiler based optimizations

if I_Tech['REAL'] == 'float':
    I_Tech['optimizations'] = ' -cl-mad-enable -cl-no-signed-zeros -cl-finite-math-only -cl-single-precision-constant'
else:
    I_Tech['optimizations'] = ' -cl-mad-enable -cl-no-signed-zeros -cl-finite-math-only'

I_RunOps['periodicx'] = 'None' #'USE_PERIODIC_X'
I_RunOps['periodicy'] = 'None' #'USE_PERIODIC_Y' 
I_RunOps['periodicz'] = 'None' #'USE_PERIODIC_Z'

I_RunOps['order'] = 4
I_RunOps['conservation_laws'] = 'ideal_MHD'
I_RunOps['testcase'] = 'alfven_periodic'
I_RunOps['plot_numerical_solution'] = 'z'
I_RunOps['save_integrals_over_time'] = False
I_RunOps['norm'] = 'LInf'

## Initialize

field_u1, field_u2 = initialize()


print('Testcase: ' + I_RunOps['testcase'] + '\nOrder: ' + str(I_RunOps['order']) + '\nTime integrator: ' + \
        I_TI['time_integrator'] +' \nDT: ' + str(I_TI['DT']) + '\nN_STEPS: ' +str(I_TI['num_steps']) + " FINAL_TIME: " + \
        str(I_TI['final_time']) + '\nDX: ' + str(I_Mesh['DX']) + ' NODES_X: ' +str(I_Mesh['NODES_X'])+'\nDY: '+ str(I_Mesh['DY']) + \
        ' NODES_Y: ' + str(I_Mesh['NODES_Y']) + '\nDZ: ' + str(I_Mesh['DZ'])+  ' NODES_Z: ' + str(I_Mesh['NODES_Z']) + '\nREAL: ' + \
        str(I_Tech['REAL']))

compute_numerical_solution(field_u1, field_u2)

rel_err = I_Results['rel_err']
for comp in range(I_BalanceLaws['NUM_CONSERVED_VARS']):
    print("Relative Error of Field Component " + str(comp) +": " + str(100*rel_err[comp]))


field_u1_reshaped = np.reshape(field_u1, (I_Tech['NUM_NODES_PAD'], I_BalanceLaws['NUM_TOTAL_VARS']))


save_all_variables(field_u1_reshaped, "output.hdf5", ['rho', 'px', 'py', 'pz', 'E', 'Bx', 'By', 'Bz'])
plot_2D(field_u1_reshaped, I_RunOps['plot_numerical_solution'], I_Mesh['NODES_X'], I_Mesh['NODES_Y'], I_Mesh['NODES_Z'], 'Numerical Solution', 1,1)
