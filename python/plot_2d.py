import matplotlib.pyplot as plt
import math

def plot_2D(u, CSs, NODES_X, NODES_Y, NODES_Z, plot_title, field_idx1, field_idx2):
    numNodes = NODES_X*NODES_Y*NODES_Z
    field = u[:,field_idx1-1:field_idx2] #Matlab counts from 1
    mag = np.zeros((numNodes, 1))
    if field.shape[1] == 1:
        mag[:-1] = field[1:numNodes]
    else:
        for i in range(numNodes):
            mag[i] = np.linalg.norm(field[i, :])
            print("Warning. This is only the absolute value")

    field_cs = np.reshape(mag, (NODES_Z, NODES_Y, NODES_X))

    for CS in CSs:
        if CS[0] == 'x':
            title = plot_title + ' X-Cross-Section'
            v = field_cs[:,:,int(math.floor(NODES_X)/2)]
            xtext = 'Y'
            ytext = 'Z'
        if CS[0] == 'y':
            title = plot_title + ' Y-Cross-Section'
            v = field_cs[:, int(math.floor(NODES_Y)/2)]
            xtext = 'X'
            ytext = 'Z'
        if CS[0] =='z':
            ytext = 'Y'
            xtext = 'X'
            title = plot_title + ' Z-Cross-Section'
            v = field_cs[int(math.floor(NODES_Z)/2),:,:]

        a = plt.pcolor(v)
        plt.colorbar(a)
        plt.title(title)
        plt.xlabel(xtext + " Axis (Nodes)")
        plt.ylabel(ytext + " Axis (Nodes)")
    plt.show()

def plot_fieldlines_B(u, CSs, NODES_X, NODES_Y, NODES_Z, plot_title, field_idx1, field_idx2):
    CS = CSs[0]
    U = u[:, field_idx1-1]
    V = u[:, field_idx2-1]
    U = np.reshape(U, (NODES_Z, NODES_Y, NODES_X))
    V = np.reshape(V, (NODES_Z, NODES_Y, NODES_X))
    if CS[0] == 'x':
        a = NODES_Y
        b = NODES_Z
        U = U[:,:,int(math.floor(NODES_X)/2)]
        V = V[:, :, int(math.floor(NODES_X)/2)]

    if CS[0] == 'y':
        a = NODES_X
        b = NODES_Z
        U = U[:, int(math.floor(NODES_Y)/2),:]
        V = V[:, int(math.floor(NODES_Y)/2),:]
    if CS[0] == 'z':
        a = NODES_X
        b = NODES_Y
        U = U[int(math.floor(NODES_Z)/2), :, :]
        V = V[int(math.floor(NODES_Z)/2), :, :]

    Y, X = np.mgrid[0:a, 0:b]
    plt.streamplot(X, Y, U, V)
    plt.show()
