import h5py

def save_all_variables(array, filename, dataset_names):
    print(array.shape)
    NODES_X = I_Mesh['NODES_X']
    NODES_Y = I_Mesh['NODES_Y']
    NODES_Z = I_Mesh['NODES_Z']
    #array = np.reshape(field_u1, (I_Tech['NUM_NODES_PAD'], I_BalanceLaws['NUM_TOTAL_VARS']))
    f = h5py.File(filename, "w")
    for i in range(len(dataset_names)):
        con_var = np.reshape(array[0:NODES_X*NODES_Y*NODES_Z, i], (NODES_Z, NODES_Y, NODES_X))
        f.create_dataset(dataset_names[i], data = con_var)

    for key in I_Tech:
        f.attrs[key] = I_Tech[key]
    for key in I_Mesh:
        f.attrs[key] = I_Mesh[key]
    for key in I_TI:
        f.attrs[key] = I_TI[key]
    for key in I_RunOps:
        f.attrs[key] = I_RunOps[key]
    for key in I_BalanceLaws:
        f.attrs[key] = I_BalanceLaws[key]


         
