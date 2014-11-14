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

    # identify cells sharing same edges as a cell
    jcell = np.zeros((nj, 3))
    for j in range(0, nj):
        k1 = jedge[j, 0]
        k2 = jedge[j, 1]
        k3 = jedge[j, 2]
        jcell[j,0] = kcell[k1,1] if (kcell[k1,0] == j) else kcell[k1,0]
        jcell[j,1] = kcell[k2,1] if (kcell[k2,0] == j) else kcell[k2,0]
        jcell[j,2] = kcell[k3,1] if (kcell[k3,0] == j) else kcell[k3,0]

    # compute cell areas
    area = np.zeros(nj)
    for j in range(0, nj):
        node1 = jnode[j, 0]
        node2 = jnode[j, 1]
        node3 = jnode[j, 2]
        x1 = xn[node1, 0]; y1 = xn[node1, 1]
        x2 = xn[node2, 0]; y2 = xn[node2, 1]
        x3 = xn[node3, 0]; y3 = xn[node3, 1]

        area[j] = .5 * abs(x1*y2 + x2*y3 + x3*y1 \
                         - x2*y1 - x3*y2 - x1*y3)

    return (xn, itype, jnode, knode, ktype, kcell, jedge, gcell, jcell, area)

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
# compute speed of sound
def sound(p, rho):
    if (p*rho) <= 0:
        return np.sqrt(abs(gamma*p/rho))
        #return cmin
    else:
        return np.sqrt(abs(gamma*p/rho))

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
# evulation of the dissipation vectors
def dissp():
    # initialize dissipation
    flux_av = np.zeros((lmax, nj))

    kappa2 = visc/4.
    kappa4 = visc/256.

    # compute flow conditions at cell-centers
    pr = np.zeros(nj)
    # interior cells
    for j in range(0, nj):
        rho = q[0,j]
        u   = q[1,j] / q[0,j]
        v   = q[2,j] / q[0,j]
        et  = q[3,j]
        pr[j] = pressure(rho, u, v, et)
    # ghost cells
    pr_ghost = np.zeros(ng)
    for g in range(0, ng):
        rho = q_ghost[0,g]
        u   = q_ghost[1,g] / q_ghost[0,g]
        v   = q_ghost[2,g] / q_ghost[0,g]
        et  = q_ghost[3,g]
        pr_ghost[g] = pressure(rho, u, v, et)

    # compute undivided pseudo-Laplacians at cell-centers
    dp = np.zeros(nj)
    dq = np.zeros((lmax, nj))
    for j in range(0, nj):
        # define cells sharing same edge as cell j
        j1 = jcell[j, 0]
        j2 = jcell[j, 1]
        j3 = jcell[j, 2]
        # compute undivided pseudo-Laplaican of pressure
        p1 = pr[j1] if (j1 >= 0) else pr_ghost[-j1-1]
        p2 = pr[j2] if (j2 >= 0) else pr_ghost[-j2-1]
        p3 = pr[j3] if (j3 >= 0) else pr_ghost[-j3-1]

        dpnum = (p1-pr[j])   + (p2-pr[j])   + (p3-pr[j])
        dpden = (p1+pr[j])/2 + (p2+pr[j])/2 + (p3+pr[j])/2
        dp[j] = dpnum / dpden

        # compute undivided pseudo-Laplacian of U
        q1 = q[:,j1] if (j1 >= 0) else q_ghost[:,-j1-1]
        q2 = q[:,j2] if (j2 >= 0) else q_ghost[:,-j2-1]
        q3 = q[:,j3] if (j3 >= 0) else q_ghost[:,-j3-1]
        dq[:,j] = (q1-q[:,j]) + (q2-q[:,j]) + (q3-q[:,j])

    # collect artificial dissiplation flux
    for k in range(0, nk):
        # two cells sharing edge k
        j1 = kcell[k,0]
        j2 = kcell[k,1]
        # node points forming edge k
        i1 = knode[k,0]
        i2 = knode[k,1]
        # compute vector length
        dx1 = xn[i2, 0] - xn[i1, 0]
        dx2 = xn[i2, 1] - xn[i1, 1]
        # compute values along edge
        if (ktype[k] == 0):
            rho = .5 * (q[0,j1] + q[0,j2])
            u   = .5 * (q[1,j1] / q[0,j1] + \
                        q[1,j2] / q[0,j2])
            v   = .5 * (q[2,j1] / q[0,j1] + \
                        q[2,j2] / q[0,j2])
            et  = .5 * (q[3,j1] + q[3,j2])

            p    = pressure(rho, u, v, et)
            c    = sound(p, rho)
            lam  = abs((u*dx2 - v*dx1) / np.sqrt(dx2**2 + dx1**2)) + c
            eps2 = kappa2 * (dp[j1] + dp[j2])
            eps4 = max(0., kappa4 - eps2)
            dx   = .5 * np.sqrt(area[j1] + area[j2])

            # incorrect sign in note
            flux1 = (eps2*( q[:,j2] -  q[:,j1]) - \
                    eps4*(dq[:,j2] - dq[:,j1])) * lam * dx
        else:
            # boundary edge (no AV along boundaries since
            # mass should not be introudued from boundaries)
            flux1 = np.zeros(lmax)

        # collect flux for left and right cells
        flux_av[:,j1] += flux1
        flux_av[:,j2] -= flux1
    return (flux_av, pr)

# -----------------------------------------------------------------------------
# evulation of source term
def source():
    # defined at cell-centers
    h = np.zeros((lmax, nj))
    for j in range(0, nj):
        # node points forming cell j
        node1 = jnode[j, 0]
        node2 = jnode[j, 1]
        node3 = jnode[j, 2]
        # compute flow variables
        yloc  = (xn[node1,1] + xn[node2,1] + xn[node3,1]) / 3.
        rho   = q[0,j]
        u     = q[1,j] / q[0,j]
        v     = q[2,j] / q[0,j]
        et    = q[3,j]
        p     = pressure(rho, u, v, et)
        # compute source
        h[0,j] = - rho*v / yloc
        h[1,j] = - rho*u*v / yloc
        h[2,j] = - rho*v*v / yloc
        h[3,j] = - (et+p)*v / yloc
    return h

# -----------------------------------------------------------------------------
# evulation of residual vectors
def residual():
    res = np.zeros((lmax, nj))
    for j in range(0, nj):
        # compute residual vector
        res[:,j] = - (flux_phys[:,j] - flux_av[:,j])/area[j] + h[:,j]
    return res

# -----------------------------------------------------------------------------
# compute time step
# coarse method, need to improve
def step():
    dt = np.zeros(nj)
    for j in range(0, nj):
        dx  = .5 * np.sqrt(2*area[j])
        rho = q[0,j]
        u   = q[1,j] / q[0,j]
        v   = q[2,j] / q[0,j]
        et  = q[3,j]
        p   = pressure(rho, u, v, et)
        c   = sound(p, rho)

        dt[j] = cfl / ((abs(u) + abs(v))/dx + c*np.sqrt(2/dx**2))
    return dt

# -----------------------------------------------------------------------------
# update q in ghost cells along solid boundary
def bc_wall():
    global q_ghost
    for g in range(0,ng):
        # boundary edge
        k = gcell[g]
        # normal cell sharing edge k with the ghost cell
        j = max(kcell[k,0], kcell[k,1])

        if (ktype[k] == 1):
            # node points forming edge k
            i1 = knode[k,0]
            i2 = knode[k,1]
            # compute vector length
            dx1 = xn[i2, 0] - xn[i1, 0]
            dx2 = xn[i2, 1] - xn[i1, 1]
            dx  = np.sqrt(dx1**2 + dx2**2)

            u         = q[1,j] / q[0,j]
            v         = q[2,j] / q[0,j]
            vn        = (u*dx2 - v*dx1) / dx
            u_ghost   = u - 2 * vn * dx2 / dx
            v_ghost   = v + 2 * vn * dx1 / dx
            rho_ghost = q[0,j]
            p_ghost   = pr[j]
            et_ghost  = p_ghost/(gamma-1) + .5*rho_ghost* \
                        (u_ghost**2 + v_ghost**2)

            q_ghost[0,g] = rho_ghost
            q_ghost[1,g] = rho_ghost * u_ghost
            q_ghost[2,g] = rho_ghost * v_ghost
            q_ghost[3,g] = et_ghost
    return

# -----------------------------------------------------------------------------
# update q in ghost cells along symmetric boundary
def bc_symmetric():
    global q_ghost

    for g in range(0,ng):
        # boundary edge
        k = gcell[g]
        # normal cell sharing edge k with the ghost cell
        j = max(kcell[k,0], kcell[k,1])

        if (ktype[k] == 2):
            q_ghost[0,g] =   q[0,j]
            q_ghost[1,g] =   q[1,j]
            q_ghost[2,g] =  -q[2,j]
            q_ghost[3,g] =   q[3,j]
    return

# -----------------------------------------------------------------------------
# update q in ghost cells along inlet boundary
def bc_inflow():
    global q_ghost

    for g in range(0,ng):
        # boundary edge
        k = gcell[g]
        # normal cell sharing edge k with the ghost cell
        j = max(kcell[k,0], kcell[k,1])

        # supersonic inflow (ibcin = 2)
        if (ibcin == 2):
            m = 1.5
            if (ktype[k] == 3):
                term = 1 / (1 + .5*(gamma-1) * m**2)
                p    = term**(gamma/(gamma-1))
                rho  = p/term
                c    = sound(p, rho)
                u    = m * c
                v    = 0.
                et   = p/(gamma-1) + .5*rho*(u**2 + v**2)

                q_ghost[0,g] = rho
                q_ghost[1,g] = rho*u
                q_ghost[2,g] = rho*v
                q_ghost[3,g] = et
    return

# -----------------------------------------------------------------------------
# update q in ghost cells along outlet boundary
def bc_outflow():
    global q_ghost

    for g in range(0,ng):
        # boundary edge
        k = gcell[g]
        # normal cell sharing edge k with the ghost cell
        j = max(kcell[k,0], kcell[k,1])

        # supersonic inflow (ibcin = 2)
        if (ibcout == 2):
            if (ktype[k] == 4):
                q_ghost[:,g] = q[:,j]
    return
# -----------------------------------------------------------------------------
def calc_mach():
    global line1, line2
    sym_mach = []
    sol_mach = []
    sym_x    = []
    sol_x    = []
    for g in range(0,ng):
        # boundary edge
        k = gcell[g]
        # normal cell sharing edge k with the ghost cell
        j = max(kcell[k,0], kcell[k,1])
        
        if (ktype[k] == 1):
            rho = q[0,j]
            u   = q[1,j] / q[0,j]
            v   = q[2,j] / q[0,j]
            et  = q[3,j]
            p   = pressure(rho, u, v, et)
            c   = sound(p, rho)
            sol_mach.append(u/c)

            node1 = jnode[j, 0]
            node2 = jnode[j, 1]
            node3 = jnode[j, 2]
            xloc  = (xn[node1,0] + xn[node2,0] + xn[node3,0]) / 3.
            sol_x.append(xloc)

        if (ktype[k] == 2):
            rho = q[0,j]
            u   = q[1,j] / q[0,j]
            v   = q[2,j] / q[0,j]
            et  = q[3,j]
            p   = pressure(rho, u, v, et)
            c   = sound(p, rho)
            sym_mach.append(u/c)
        
            node1 = jnode[j, 0]
            node2 = jnode[j, 1]
            node3 = jnode[j, 2]
            xloc  = (xn[node1,0] + xn[node2,0] + xn[node3,0]) / 3.
            sym_x.append(xloc)
    if (t == 0):
        line1, = plt.plot(sym_x, sym_mach, 'o-')
        line2, = plt.plot(sol_x, sol_mach, 'o-')
        ax.set_ylim(1,3)
        fig.show()
    else:
        line1.set_ydata(sym_mach)
        line2.set_ydata(sol_mach)
        fig.canvas.draw()

    return

# -----------------------------------------------------------------------------
# four-stage Runge-Kutta scheme
def rk4():
    global q, q_ghost, flux_phys, flux_av, pr, h, res

    rk    = [1./4, 1./3, 1./2, 1.]
    q_old = np.copy(q)

    for m in range(0, 4):
        # compute flux
        flux_phys = flux()
        # compute aritificial viscosity flux
        (flux_av, pr) = dissp()
        # compute source
        h = source()
        # compute residual
        res = residual()
        # compute time step
        #dt = step()
        dt = np.ones(nj)*.005

        # update solution
        for j in range(0, nj):
            q[:,j] = q_old[:,j] + rk[m]*dt[j]*res[:,j]

        # boundary conditions
        bc_wall()
        bc_symmetric()
        bc_inflow()
        bc_outflow()

    calc_mach()
    return

# -----------------------------------------------------------------------------
# main program

# grid points (structured mesh)
mx = 65; my = 17 
# mx = 3; my = 3

# parameters
lmax = 4
gamma = 1.4
alpha = 1       # axisymmetric flow

cfl  = .1
visc = 4
eps  = 1e-5
cmin = 1e-2
tmax = 100000

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
filename = 'ft03.dat'
(xs, ys) = read_mesh(filename)

(xn, itype, jnode, knode, ktype, kcell, jedge, gcell, jcell, area) = transform()
# plot_mesh()

(q, q_ghost) = ic()

fig = plt.figure()
ax  = fig.add_subplot(111, aspect='equal')

for t in range(0,tmax):
    stop = False
    rk4()
    print 'step ' + str(t), np.amax(abs(res[0,:]))
    if (stop): break
