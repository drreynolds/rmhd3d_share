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
def load_prec():
    """Returns hprec, vprec, combo, nwaves, order, ADIcorr from prec.inp"""
    import shlex
    # load the preconditioner options
    f = open("prec.inp")
    hprec=0;  vprec=0;  method=0;  NWaves=0;  ordering=0;  ADIcorr=0;
    for line in f:
        txt = shlex.split(line)
        if ("hprec" in txt):
            hprec = int(txt[len(txt)-1].replace(',',''))
        elif ("vprec" in txt):
            vprec = int(txt[len(txt)-1].replace(',',''))
        elif ("method" in txt):
            method = int(txt[len(txt)-1].replace(',',''))
        elif ("NWaves" in txt):
            NWaves = int(txt[len(txt)-1].replace(',',''))
        elif ("ordering" in txt):
            ordering = int(txt[len(txt)-1].replace(',',''))
        elif ("ADIcorrection" in txt):
            ADIcorr = int(txt[len(txt)-1].replace(',',''))
    f.close()
    return [hprec, vprec, method, NWaves, ordering, ADIcorr]

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
def load_vals(tdump):
    """Returns x,y,z,rho,u,v,w,bx,by,bz,p,dB,j,te from a given data dump"""
    import shlex
    import numpy as np
    sdump = repr(tdump).zfill(6)
    outfile = 'output.001.' + sdump
    nx,ny,nz,xprocs,yprocs,zprocs,xbc,ybc,zbc = load_mesh()
    f = open(outfile)
    x  = np.zeros(nx*ny*nz,dtype=float); y   = np.zeros(nx*ny*nz,dtype=float); 
    z  = np.zeros(nx*ny*nz,dtype=float); rho = np.zeros(nx*ny*nz,dtype=float); 
    u  = np.zeros(nx*ny*nz,dtype=float); v   = np.zeros(nx*ny*nz,dtype=float); 
    w  = np.zeros(nx*ny*nz,dtype=float); bx  = np.zeros(nx*ny*nz,dtype=float); 
    by = np.zeros(nx*ny*nz,dtype=float); bz  = np.zeros(nx*ny*nz,dtype=float); 
    p  = np.zeros(nx*ny*nz,dtype=float); dB  = np.zeros(nx*ny*nz,dtype=float); 
    j  = np.zeros(nx*ny*nz,dtype=float); te  = np.zeros(nx*ny*nz,dtype=float); 
    idx = 0;
    for line in f:
        txt = shlex.split(line)
        if (("#" not in txt) and (txt != [])):
            x[idx]   = float(txt[0])
            y[idx]   = float(txt[1])
            z[idx]   = float(txt[2])
            rho[idx] = float(txt[3])
            u[idx]   = float(txt[4])
            v[idx]   = float(txt[5])
            w[idx]   = float(txt[6])
            bx[idx]  = float(txt[7])
            by[idx]  = float(txt[8])
            bz[idx]  = float(txt[9])
            p[idx]   = float(txt[10])
            dB[idx]  = float(txt[11])
            j[idx]   = float(txt[12])
            te[idx]  = float(txt[13])
            idx += 1
    f.close()
    # reshape 1D arrays to 3D
    xarr = np.array(x[0:nx]);
    yarr = np.zeros(ny);
    zarr = np.zeros(nz);
    idx = 0;
    for k in range(ny):
        yarr[idx] = y[k*nx];
        idx += 1;
    idx = 0;
    for k in range(nz):
        zarr[idx] = z[k*nx*ny];
        idx += 1;
    rho3D = np.reshape(rho, (nx,ny,nz), order='F')
    u3D   = np.reshape(u,   (nx,ny,nz), order='F')
    v3D   = np.reshape(v,   (nx,ny,nz), order='F')
    w3D   = np.reshape(w,   (nx,ny,nz), order='F')
    bx3D  = np.reshape(bx,  (nx,ny,nz), order='F')
    by3D  = np.reshape(by,  (nx,ny,nz), order='F')
    bz3D  = np.reshape(bz,  (nx,ny,nz), order='F')
    p3D   = np.reshape(p,   (nx,ny,nz), order='F')
    dB3D  = np.reshape(dB,  (nx,ny,nz), order='F')
    j3D   = np.reshape(j,   (nx,ny,nz), order='F')
    te3D  = np.reshape(te,  (nx,ny,nz), order='F')
    return [xarr,yarr,zarr,rho3D,u3D,v3D,w3D,bx3D,by3D,bz3D,p3D,dB3D,j3D,te3D]
##########

