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

    # edge type:
    # 0 - interior; 1 - solid boundary; 2 - symmetric line
    # 3 - inflow;   4 - outflow
    ktype = np.zeros(nk)

    # symmetric line
    j = 0
    for i in range(0, mx-1):
        k = mx * j + i
        edge1 = [k, k+1]
        edge2 = [k+1, k]
        if (edge2 in klist): edge1 = edge2
        ktype[klist.index(edge1)] = 2
    # solid boundary
    j = my-1
    for i in range(0, mx-1):
        k = mx * j + i
        edge1 = [k, k+1]
        edge2 = [k+1, k]
        if (edge2 in klist): edge1 = edge2
        ktype[klist.index(edge1)] = 1
    # inflow
    i = 0
    for j in range(0, my-1):
        k = mx * j + i
        edge1 = [k, k+mx]
        edge2 = [k+mx, k]
        if (edge2 in klist): edge1 = edge2
        ktype[klist.index(edge1)] = 3
    # outflow
    i = mx-1
    for j in range(0, my-1):
        k = mx * j + i
        edge1 = [k, k+mx]
        edge2 = [k+mx, k]
        if (edge2 in klist): edge1 = edge2
        ktype[klist.index(edge1)] = 4

    # identify cells sharing edge k
    # -1 means no cell (ghost cell) is on the side
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

    # other cell data
    # identify 3 edges forming the cells
    jedge = np.zeros((nj, 3))
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
        
        jedge[j, 0] = k1
        jedge[j, 1] = k2
        jedge[j, 2] = k3

    # ghost cells (defined by the sharing edge with normal cells)
    gcelllist = []
    for k in range(0, nk):
        if (ktype[k] != 0):
            gcelllist.append(k)
    gcell = np.array(gcelllist)

    # negative number in kcell identify corresponding ghost cell
    for g in range(0, ng):
        if kcell[gcell[g], 0] == -1:
            kcell[gcell[g], 0] = -g-1
        else:
            kcell[gcell[g], 1] = -g-1

    return (xn, itype, jnode, knode, ktype, kcell, jedge, gcell)

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
# compute pressure
def pressure(rho, u, v, et):
    return (gamma-1) * (et - .5 * rho * (u**2 + v**2))

# -----------------------------------------------------------------------------
def plot_mesh():
    fig = plt.figure()
    ax  = fig.add_subplot(111, aspect='equal') 
    for j in range(0, nj):
        cycle = np.append(jnode[j,:], jnode[j,0]).tolist()
        ax.plot(xn[cycle, 0], xn[cycle,1], '-o', color='k')
    plt.show()
    return
# -----------------------------------------------------------------------------
# initialize conservative variables
def ic():
    # defined at cell-centers
    q = np.ones((lmax, nj))
    if (ibcin == 1):
        q[0,:] = 1.
        q[1,:] = .5
        q[2,:] = 0.
        q[3,:] = 1.
    if (ibcin == 2):
        q[0,:] = 1.
        q[1,:] = 1.5
        q[2,:] = 0.
        q[3,:] = 2.

    # ghost cells
    q_ghost = np.ones((lmax, ng))
    if (ibcin == 1):
        q_ghost[0,:] = 1.
        q_ghost[1,:] = .5
        q_ghost[2,:] = 0.
        q_ghost[3,:] = 1.
    if (ibcin == 2):
        q_ghost[0,:] = 1.
        q_ghost[1,:] = 1.5
        # q_ghost[1,:] = 1.
        q_ghost[2,:] = 0.
        q_ghost[3,:] = 2.
    return q, q_ghost

# -----------------------------------------------------------------------------
# evulation of flux
def flux():
    # initialize flux
    flux_phys = np.zeros((lmax, nj))

    for k in range(0, nk):
        # node points forming edge k
        i1 = knode[k,0]
        i2 = knode[k,1]
        # compute vector length
        dx1 = xn[i2, 0] - xn[i1, 0]
        dx2 = xn[i2, 1] - xn[i1, 1]
        # two cells sharing edge k
        j1 = kcell[k,0]
        j2 = kcell[k,1]
        # compute average flow quantities along edge
        if (ktype[k] == 0):
            rho = .5 * (q[0,j1] + q[0,j2])
            u   = .5 * (q[1,j1] / q[0,j1] + \
                        q[1,j2] / q[0,j2])
            v   = .5 * (q[2,j1] / q[0,j1] + \
                        q[2,j2] / q[0,j2])
            et  = .5 * (q[3,j1] + q[3,j2])
            p   = pressure(rho, u, v, et)
        elif (j1 < 0):
            rho = .5 * (q_ghost[0,-j1-1] + q[0,j2])
            u   = .5 * (q_ghost[1,-j1-1] / q_ghost[0,-j1-1] + \
                        q[1,j2] / q[0,j2])
            v   = .5 * (q_ghost[2,-j1-1] / q_ghost[0,-j1-1] + \
                        q[2,j2] / q[0,j2])
            et  = .5 * (q_ghost[3,-j1-1] + q[3,j2])
            p   = pressure(rho, u, v, et)
        elif (j2 < 0):
            rho = .5 * (q_ghost[0,-j2-1] + q[0,j1])
            u   = .5 * (q_ghost[1,-j2-1] / q_ghost[0,-j2-1] + \
                        q[1,j1] / q[0,j1])
            v   = .5 * (q_ghost[2,-j2-1] / q_ghost[0,-j2-1] + \
                        q[2,j1] / q[0,j1])
            et  = .5 * (q_ghost[3,-j2-1] + q[3,j1])
            p   = pressure(rho, u, v, et)
        # compute flux vector along edge
        flux1 = rho*(u*dx2 - v*dx1)
        flux2 = (rho*u**2 + p)*dx2 - (rho*u*v)*dx1
        flux3 = (rho*u*v)*dx2 - (rho*v**2 + p)*dx1
        flux4 = (et+p)*(u*dx2 - v*dx1)
        # collect flux for left and right cells
        if (ktype[k] == 0):
            flux_phys[0,j1] += flux1
            flux_phys[0,j2] -= flux1
            flux_phys[1,j1] += flux2
            flux_phys[1,j2] -= flux2
            flux_phys[2,j1] += flux3
            flux_phys[2,j2] -= flux3
            flux_phys[3,j1] += flux4
            flux_phys[3,j2] -= flux4
        elif (j1 < 0):
            flux_phys[0,j2] -= flux1
            flux_phys[1,j2] -= flux2
            flux_phys[2,j2] -= flux3
            flux_phys[3,j2] -= flux4
        elif (j2 < 0):
            flux_phys[0,j1] += flux1
            flux_phys[1,j1] += flux2
            flux_phys[2,j1] += flux3
            flux_phys[3,j1] += flux4
    return flux_phys

# -----------------------------------------------------------------------------
# main program

# grid points (structured mesh)
# mx = 65; my = 17 
mx = 3; my = 3

# parameters
lmax = 4
gamma = 1.4
alpha = 1       # axisymmetric flow

# ni - number of points; nj - number of cells;
# nk - number of edges;  np - number of ghost cells
ni = mx * my
nj = (mx-1) * (my-1) * 2
nk = (mx-1) * my + (my-1) * mx + (mx-1) * (my-1)
ng = (mx + my - 2) * 2

# boundary conditions
ibcin  = 2
ibcout = 2

# coordinates of structured mesh
filename = 'test.dat'
(xs, ys) = read_mesh(filename)

(xn, itype, jnode, knode, ktype, kcell, jedge, gcell) = transform()
# plot_mesh()

(q, q_ghost) = ic()

flux_phys = flux()
print flux_phys
