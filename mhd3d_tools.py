# matplotlib-based plotting utilities for RMHD3d tests
# Daniel R. Reynolds, reynolds@smu.edu


##########
def load_cons():
    """Returns times, energy, mass arrays from the conservation.txt file"""
    import shlex
    import numpy as np
    # load the conservation.txt file
    # (columns are:  time, total energy, total mass)
    f = open("conservation.txt")
    times=[];  tenergy=[];  totmass=[]; 
    for line in f:
        txt = shlex.split(line)
        times.append(float(txt[0]));
        tenergy.append(float(txt[1]));
        totmass.append(float(txt[2]));
    f.close()
    return [np.array(times), np.array(tenergy), np.array(totmass)]

##########
def load_energies():
    """Returns kenergy, menergy, ienergy from the energyHistory.txt file"""
    import shlex
    import numpy as np
    # load the energyHistory.txt file
    # (columns are:  time, kinetic energy, magnetic energy, internal energy)
    f = open("energyHistory.txt")
    kenergy=[];  menergy=[];  ienergy=[]; 
    for line in f:
        txt = shlex.split(line)
        kenergy.append(float(txt[1]));
        menergy.append(float(txt[2]));
        ienergy.append(float(txt[3]));
    f.close()
    return [np.array(kenergy), np.array(menergy), np.array(ienergy)]

##########
def load_props():
    """Returns mu, eta, kappa, gamma from prop.inp"""
    import shlex
    # load the physical simulation properties
    f = open("prop.inp")
    mu=0.0;  eta=0.0;  kappa=0.0;  gamma=5.0/3.0;
    for line in f:
        txt = shlex.split(line)
        if ("mu" in txt):
            mu = float(txt[len(txt)-1].replace('d','e').replace(',',''))
        elif ("eta" in txt):
            eta = float(txt[len(txt)-1].replace('d','e').replace(',',''))
        elif ("kappa" in txt):
            kappa = float(txt[len(txt)-1].replace('d','e').replace(',',''))
        elif ("gamma" in txt):
            gamma = float(txt[len(txt)-1].replace('d','e').replace(',',''))
    f.close()
    return [mu, eta, kappa, gamma]    

##########
def load_mhd():
    """Returns xl,xr,yl,yr,zl,zr,ndump from mhd.inp"""
    import shlex
    # load the grid bounds
    f = open("mhd.inp")
    xl=0.0;  xr=1.0;  yl=0.0;  yr=1.0;  zl=0.0;  zr=1.0;  ndump=1;
    for line in f:
        txt = shlex.split(line)
        if ("xl" in txt):
            xl = float(txt[len(txt)-1].replace('d','e').replace(',',''))
        elif ("xr" in txt):
            xr = float(txt[len(txt)-1].replace('d','e').replace(',',''))
        elif ("yl" in txt):
            yl = float(txt[len(txt)-1].replace('d','e').replace(',',''))
        elif ("yr" in txt):
            yr = float(txt[len(txt)-1].replace('d','e').replace(',',''))
        elif ("zl" in txt):
            zl = float(txt[len(txt)-1].replace('d','e').replace(',',''))
        elif ("zr" in txt):
            zr = float(txt[len(txt)-1].replace('d','e').replace(',',''))
        elif ("ndump" in txt):
            ndump = int(txt[len(txt)-1].replace(',',''))
    f.close()
    return [xl, xr, yl, yr, zl, zr, ndump]

##########
def load_mesh():
    """Returns nx,ny,nz,xprocs,yprocs,zprocs,xbc,ybc,zbc from mesh.inp"""
    import shlex
    # load the mesh information
    f = open("mesh.inp")
    nx=1; ny=1; nz=1; xbc=0; ybc=0; zbc=0;
    xprocs=1; yprocs=1; zprocs=1; 
    for line in f:
        txt = shlex.split(line)
        if ("nx" in txt):
            nx = int(txt[len(txt)-1].replace(',',''))
        elif ("ny" in txt):
            ny = int(txt[len(txt)-1].replace(',',''))
        elif ("nz" in txt):
            nz = int(txt[len(txt)-1].replace(',',''))
        elif ("xprocs" in txt):
            xprocs = int(txt[len(txt)-1].replace(',',''))
        elif ("yprocs" in txt):
            yprocs = int(txt[len(txt)-1].replace(',',''))
        elif ("zprocs" in txt):
            zprocs = int(txt[len(txt)-1].replace(',',''))
        elif ("xbc" in txt):
            xbc = int(txt[len(txt)-1].replace(',',''))
        elif ("ybc" in txt):
            ybc = int(txt[len(txt)-1].replace(',',''))
        elif ("zbc" in txt):
            zbc = int(txt[len(txt)-1].replace(',',''))
    f.close()
    return [nx, ny, nz, xprocs, yprocs, zprocs, xbc, ybc, zbc];

##########
def load_RunHistory():
    """Returns Newton, Krylov, Fnorm, divB, SolTime from Run.history"""
    import shlex
    import numpy as np
    # load the Run history information
    f = open("Run.history")
    Newton=[];  Krylov=[];  Fnorm=[];  divB=[];  SolTime=0.0;
    for line in f:
        txt = shlex.split(line)
        if ("Num. Newton Iterations" in txt):
            Newton.append(int(txt[len(txt)-1]))
        elif ("Number Krylov Iterations" in txt):
            Krylov.append(int(txt[len(txt)-1]))
        elif ("Final scaled norm of f" in txt):
            Fnorm.append(float(txt[len(txt)-1]).replace('d','e'))
        elif ("divB" in txt):
            divB = float(txt[len(txt)-4].replace('d','e'))
        elif ("Total Solution Time" in txt):
            SolTime = float(txt[len(txt)-1].replace('d','e'))
    f.close()
    return [np.array(Newton), np.array(Krylov), np.array(Fnorm), 
            np.array(divB), np.array(SolTime)];
    
##########
def load_proclayout():
    """Returns nprocs and layout, where layout is a 2D array of processor layout"""
    """(each row is 1 proc, columns are processor location in x,y,z directions)"""
    """from proc_layout.txt"""
    import shlex
    import numpy as np
    # load the mesh information
    data = np.loadtxt("proc_layout.txt")
    xprocs = int(data[0,1])
    yprocs = int(data[0,2])
    zprocs = int(data[0,3])
    nprocs = xprocs*yprocs*zprocs
    layout = np.zeros((nprocs,3), dtype=int)
    for iproc in range(nprocs):
        layout[iproc,0] = data[iproc+1,1]
        layout[iproc,1] = data[iproc+1,2]
        layout[iproc,2] = data[iproc+1,3]
    return [nprocs, layout];

##########
def load_mygrid(xloc,yloc,zloc):
    """Returns x,y,z arrays for the processor at process location xloc,yloc,zloc"""
    import shlex
    import numpy as np
    # load the mesh information
    xfile = 'x.grid.' + repr(xloc).zfill(3)
    yfile = 'y.grid.' + repr(yloc).zfill(3)
    zfile = 'z.grid.' + repr(zloc).zfill(3)
    xgrid = np.loadtxt(xfile)
    ygrid = np.loadtxt(yfile)
    zgrid = np.loadtxt(zfile)
    return [xgrid, ygrid, zgrid];

##########
def count_outputs():
    """Returns Nout, the total number of time dumps"""
    import shlex
    import os.path
    # get ndump
    xl,xr,yl,yr,zl,zr,ndump = load_mhd()
    maxoutputs = 10000000
    Nout = 0
    for iout in range(0,maxoutputs,ndump):
        outfile = 'data' + repr(iout).zfill(6)
        if (os.path.isdir(outfile)):
            Nout += 1
        else:
            break
    return Nout;

##########
def load_vals(tdump):
    """Returns x,y,z,rho,u,v,w,bx,by,bz,p,dB,jz,ptot from a given data dump"""
    import shlex
    import numpy as np
    sdump = repr(tdump).zfill(6)
    # load mesh information
    nx,ny,nz,xprocs,yprocs,zprocs,xbc,ybc,zbc = load_mesh()
    # load domain information
    xl,xr,yl,yr,zl,zr,ndump = load_mhd()
    # load processor layout
    nprocs,layout = load_proclayout()
    # check that layout and mesh.inp match
    if (xprocs*yprocs*zprocs != nprocs):
        print 'Error: proc_layout.txt and mesh.inp disagree on total proc count!'
        return -1
    # set mesh increments
    dx = (xr-xl)/nx
    dy = (yr-yl)/ny
    dz = (zr-zl)/nz
    # set x,y,z Cartesian grid information
    x = np.linspace(xl+0.5*dx,xr-0.5*dx,nx)
    y = np.linspace(yl+0.5*dy,yr-0.5*dy,ny)
    z = np.linspace(zl+0.5*dz,zr-0.5*dz,nz)
    # set generic local problem dimensions (last proc. in a dir. may differ)
    nxlocal = nx/xprocs
    nylocal = ny/yprocs
    nzlocal = nz/zprocs
    # initialize solution output
    rho  = np.zeros((nx,ny,nz), dtype=float, order='F')
    u    = np.zeros((nx,ny,nz), dtype=float, order='F')
    v    = np.zeros((nx,ny,nz), dtype=float, order='F')
    w    = np.zeros((nx,ny,nz), dtype=float, order='F')
    bx   = np.zeros((nx,ny,nz), dtype=float, order='F')
    by   = np.zeros((nx,ny,nz), dtype=float, order='F')
    bz   = np.zeros((nx,ny,nz), dtype=float, order='F')
    p    = np.zeros((nx,ny,nz), dtype=float, order='F')
    dB   = np.zeros((nx,ny,nz), dtype=float, order='F')
    jcur = np.zeros((nx,ny,nz), dtype=float, order='F')
    ptot = np.zeros((nx,ny,nz), dtype=float, order='F')
    # iterate over the output files
    for iproc in range(nprocs):
        # determine output file for this processor
        sproc = repr(iproc).zfill(6)
        outfile = 'data' + sdump + '/output' + sproc + '.txt'
        # determine local processor extents in global grid
        istart = (layout[iproc,0]-1)*nxlocal
        iend   = min(istart+nxlocal-1, nx-1)
        jstart = (layout[iproc,1]-1)*nylocal
        jend   = min(jstart+nylocal-1, ny-1)
        kstart = (layout[iproc,2]-1)*nzlocal
        kend   = min(kstart+nzlocal-1, nz-1)
        data = np.loadtxt(outfile)
        # load local data file
        data = np.reshape(data, (nxlocal,nylocal,nzlocal,11), order='F')
        # store local data file into global outputs
        #    rho
        for k in range(kend-kstart+1):
            for j in range(jend-jstart+1):
                for i in range(iend-istart+1):
                    rho[istart+i,jstart+j,kstart+k] = data[i,j,k,0]
        #    u
        for k in range(kend-kstart+1):
            for j in range(jend-jstart+1):
                for i in range(iend-istart+1):
                    u[istart+i,jstart+j,kstart+k] = data[i,j,k,1]
        #    v
        for k in range(kend-kstart+1):
            for j in range(jend-jstart+1):
                for i in range(iend-istart+1):
                    v[istart+i,jstart+j,kstart+k] = data[i,j,k,2]
        #    w
        for k in range(kend-kstart+1):
            for j in range(jend-jstart+1):
                for i in range(iend-istart+1):
                    w[istart+i,jstart+j,kstart+k] = data[i,j,k,3]
        #    bx
        for k in range(kend-kstart+1):
            for j in range(jend-jstart+1):
                for i in range(iend-istart+1):
                    bx[istart+i,jstart+j,kstart+k] = data[i,j,k,4]
        #    by
        for k in range(kend-kstart+1):
            for j in range(jend-jstart+1):
                for i in range(iend-istart+1):
                    by[istart+i,jstart+j,kstart+k] = data[i,j,k,5]
        #    bz
        for k in range(kend-kstart+1):
            for j in range(jend-jstart+1):
                for i in range(iend-istart+1):
                    bz[istart+i,jstart+j,kstart+k] = data[i,j,k,6]
        #    p
        for k in range(kend-kstart+1):
            for j in range(jend-jstart+1):
                for i in range(iend-istart+1):
                    p[istart+i,jstart+j,kstart+k] = data[i,j,k,7]
        #    dB
        for k in range(kend-kstart+1):
            for j in range(jend-jstart+1):
                for i in range(iend-istart+1):
                    dB[istart+i,jstart+j,kstart+k] = data[i,j,k,8]
        #    jcur
        for k in range(kend-kstart+1):
            for j in range(jend-jstart+1):
                for i in range(iend-istart+1):
                    jcur[istart+i,jstart+j,kstart+k] = data[i,j,k,9]
        #    ptot
        for k in range(kend-kstart+1):
            for j in range(jend-jstart+1):
                for i in range(iend-istart+1):
                    ptot[istart+i,jstart+j,kstart+k] = data[i,j,k,10]
    return [x,y,z,rho,u,v,w,bx,by,bz,p,dB,jcur,ptot]
##########

