import matplotlib.pyplot as plt
import math

def plot_2D(u, CSs, NODES_X, NODES_Y, NODES_Z, plot_title, field_idx1, field_idx2):
    numNodes = NODES_X*NODES_Y*NODES_Z
    field = u[:,field_idx1-1:field_idx2] #Matlab counts from 1
    mag = np.zeros((numNodes, 1))
    if field.shape[0] == 1:
        mag[:] = field[1:numNodes]
    else:
        for i in range(numNodes):
            mag[i] = np.linalg.norm(field[i, :])

    field_cs = np.reshape(mag, (NODES_Z, NODES_Y, NODES_X))

    for CS in CSs:
        if CS[0] == 'x':
            title = plot_title + 'X-Cross-Section'
            v = field_cs[:,:,int(math.floor(NODES_X)/2)]
        if CS[0] == 'y':
            title = plot_title + 'Y-Cross-Section'
            v = field_cs[:, int(math.floor(NODES_Y)/2)]
        if CS[0] =='z':
            title = plot_title + 'Z-Cross-Section'
            v = field_cs[int(math.floor(NODES_Z)/2), :, :]

        plt.contourf(v)
        plt.title(title)
    plt.show()

	
