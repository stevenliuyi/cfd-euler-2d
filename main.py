import numpy as np
import matplotlib.pyplot as plt

# -----------------------------------------------------------------------------
def read_mesh(filename):
    f = open(filename)
    x = np.zeros((mx, my))
    y = np.zeros((mx, my))
    i = 0; j = 0
    for line in f:
        x[i,j] = line.split()[0]
        y[i,j] = line.split()[1]
        if (j == my-1):
            i += 1; j = 0
        else:
            j += 1
        if (i == mx): break

    f.close()
    return (x, y)

# -----------------------------------------------------------------------------
# transform structured mesh to unstructred mesh
def transform():
    # generate node points
    xn = np.zeros((ni, 2))
    # node type:
    # 0 - interior; 1 - solid boundary; 2 - symmetric line
    # 3 - inflow;   4 - outflow
    itype = np.zeros(ni)

    for i in range(0, mx):
        for j in range(0, my):
            k = mx * j + i
            xn[k, 0] = xs[i, j]
            xn[k, 1] = ys[i, j]
            if (j == 0):    itype[k] = 2
            if (j == my-1): itype[k] = 1
            if (i == 0):    itype[k] = 3
            if (i == mx-1): itype[k] = 4

    # generate the cells
    jnode = np.zeros((nj, 3))
    for i in range(0, mx-1):
        for j in range(0, my-1):
            k = (mx-1)*j + i
            jnode[2*k,   0] = mx * j     + i
            jnode[2*k,   1] = mx * j     + (i+1)
            jnode[2*k,   2] = mx * (j+1) + (i+1)
            jnode[2*k+1, 0] = mx * j     + i
            jnode[2*k+1, 1] = mx * (j+1) + (i+1)
            jnode[2*k+1, 2] = mx * (j+1) + i
    
    return (xn, itype, jnode)

# -----------------------------------------------------------------------------
# main program

# grid points (structured mesh)
# mx = 65; my = 17 
mx = 4; my = 4

# ni - number of points; nj - number of cells; nk - number of edges
ni = mx * my
nj = (mx-1) * (my-1) * 2

# coordinates of structured mesh
filename = 'test.dat'
(xs, ys) = read_mesh(filename)

(xn, itype, jnode) = transform()
