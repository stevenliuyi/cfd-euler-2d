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

    # generate the cells (3 node points of cells)
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
    
    # generate edge data
    # 2 node points of edges
    klist = []
    for j in range(0, nj):
        edge1 = [jnode[j,0], jnode[j,1]]
        edge2 = [jnode[j,1], jnode[j,0]]
        edge3 = [jnode[j,0], jnode[j,2]]
        edge4 = [jnode[j,2], jnode[j,0]]
        edge5 = [jnode[j,1], jnode[j,2]]
        edge6 = [jnode[j,2], jnode[j,1]]
        if not((edge1 in klist) or (edge2 in klist)):
            klist.append(edge1)
        if not((edge3 in klist) or (edge4 in klist)):
            klist.append(edge3)
        if not((edge5 in klist) or (edge6 in klist)):
            klist.append(edge5)
    knode = np.array(klist)

    # identify cells sharing edge k
    # -1 means no cell is on the side
    kcell = -np.ones((nk, 2))
    for j in range(0,nj):
        edge1 = [jnode[j,0], jnode[j,1]]
        edge2 = [jnode[j,1], jnode[j,0]]
        edge3 = [jnode[j,0], jnode[j,2]]
        edge4 = [jnode[j,2], jnode[j,0]]
        edge5 = [jnode[j,1], jnode[j,2]]
        edge6 = [jnode[j,2], jnode[j,1]]
        if (edge2 in klist): edge1 = edge2
        if (edge4 in klist): edge3 = edge4
        if (edge6 in klist): edge5 = edge6

        k1 = klist.index(edge1)
        k2 = klist.index(edge3)
        k3 = klist.index(edge5)
        
        kcell[k1, side(edge1, jnode[j,2], xn)] = j
        kcell[k2, side(edge3, jnode[j,1], xn)] = j
        kcell[k3, side(edge5, jnode[j,0], xn)] = j

    return (xn, itype, jnode, knode, kcell)

# -----------------------------------------------------------------------------
# determine a cell is on the left or right of a edge
def side(edge, point, xn):
    v1 = (xn[edge[1],0] - xn[edge[0],0], \
          xn[edge[1],1] - xn[edge[0],1])
    v2 = (xn[point  ,0] - xn[edge[0],0], \
          xn[point  ,1] - xn[edge[0],1])
    if (v1[0]*v2[1] - v2[0]*v1[1]) > 0:
        return 0    # left
    else:
        return 1    # right

# -----------------------------------------------------------------------------
def plot_mesh():
    fig = plt.figure()
    ax  = fig.add_subplot(111, aspect='equal') 
    for j in range(0, nj):
        cycle = np.append(jnode[j,:], jnode[j,0]).tolist()
        ax.plot(xn[cycle, 0], xn[cycle,1], '-o', color='k')
    plt.show()
# -----------------------------------------------------------------------------
# main program

# grid points (structured mesh)
# mx = 65; my = 17 
mx = 4; my = 4

# ni - number of points; nj - number of cells; nk - number of edges
ni = mx * my
nj = (mx-1) * (my-1) * 2
nk = (mx-1) * my + (my-1) * mx + (mx-1) * (my-1)

# coordinates of structured mesh
filename = 'test.dat'
(xs, ys) = read_mesh(filename)

(xn, itype, jnode, knode, kcell) = transform()
# plot_mesh()
