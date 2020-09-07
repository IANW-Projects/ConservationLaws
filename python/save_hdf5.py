import h5py

def save_all_variables(array, filename, dataset_names):
    NODES_X = I_Mesh['NODES_X']
    NODES_Y = I_Mesh['NODES_Y']
    NODES_Z = I_Mesh['NODES_Z']
    array = np.reshape(array, (I_Tech['NUM_NODES_PAD'], I_BalanceLaws['NUM_TOTAL_VARS']))
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

    f.close()

def reload_all_variables(array, filename, dataset_names):
    NODES_X = I_Mesh['NODES_X']
    NODES_Y = I_Mesh['NODES_Y']
    NODES_Z = I_Mesh['NODES_Z']
    f = h5py.File(filename, "r")
    ar_shape = array.shape
    ar = np.reshape(array, (I_Tech['NUM_NODES_PAD'], I_BalanceLaws['NUM_TOTAL_VARS']))

    for i in range(len(dataset_names)):
        dset = f[dataset_names[i]]
        ar[0:NODES_X*NODES_Y*NODES_Z, i] = np.reshape(dset, NODES_Z * NODES_Y* NODES_Z)

    array = np.reshape(ar, ar_shape)
    return array
    
