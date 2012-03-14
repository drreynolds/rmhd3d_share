# matplotlib-based plotting script for RMHD3d tests
# Daniel R. Reynolds, reynolds@smu.edu

# imports
from pylab import *
import shlex
import mhd3d_tools as mhd3d

# load the relevant files, and parameters
times, energy, mass = mhd3d.load_cons()
xl, xr, yl, yr, zl, zr, ndump = mhd3d.load_mhd()
nx, ny, nz, xprocs, yprocs, zprocs, xbc, ybc, zbc = mhd3d.load_mesh()

# set the total number of snapshots
Nt = mhd3d.count_outputs()
print 'printing ',Nt,' solution snapshots'

# slice through the center of the domain (in z-direction)
sl = int(nz/2)

# slice through offset plane for Jz current (in x-direction)
slx = int(nx/2); 

# loop over snapshots, loading values and plotting
for tstep in range(0,Nt*ndump,ndump):
    print '  step ',tstep
    
    # load solution data
    x,y,z,rho,u,v,w,bx,by,bz,p,dB,j,ptot = mhd3d.load_vals(tstep)

    # rho contour
    figure()
    s = rho[:,:,sl]
    h = imshow(s, hold=False, extent=(yl, yr, xl, xr), origin='lower')
    colorbar(h, orientation='horizontal')
    title('rho, tstep = ' + repr(tstep))
    savefig('rho_' + repr(tstep).zfill(6) + '.png')
        
    # u contour
    figure()
    s = u[:,:,sl]
    h = imshow(s, hold=False, extent=(yl, yr, xl, xr), origin='lower')
    colorbar(h, orientation='horizontal')
    title('u, tstep = ' + repr(tstep))
    savefig('u_' + repr(tstep).zfill(6) + '.png')
        
    # v contour
    figure()
    s = v[:,:,sl]
    h = imshow(s, hold=False, extent=(yl, yr, xl, xr), origin='lower')
    colorbar(h, orientation='horizontal')
    title('v, tstep = ' + repr(tstep))
    savefig('v_' + repr(tstep).zfill(6) + '.png')
        
    # w contour
    figure()
    s = w[:,:,sl]
    h = imshow(s, hold=False, extent=(yl, yr, xl, xr), origin='lower')
    colorbar(h, orientation='horizontal')
    title('w, tstep = ' + repr(tstep))
    savefig('w_' + repr(tstep).zfill(6) + '.png')
        
    # bx contour
    figure()
    s = bx[:,:,sl]
    h = imshow(s, hold=False, extent=(yl, yr, xl, xr), origin='lower')
    colorbar(h, orientation='horizontal')
    title('bx, tstep = ' + repr(tstep))
    savefig('bx_' + repr(tstep).zfill(6) + '.png')
        
    # by contour
    figure()
    s = by[:,:,sl]
    h = imshow(s, hold=False, extent=(yl, yr, xl, xr), origin='lower')
    colorbar(h, orientation='horizontal')
    title('by, tstep = ' + repr(tstep))
    savefig('by_' + repr(tstep).zfill(6) + '.png')
        
    # bz contour
    figure()
    s = bz[:,:,sl]
    h = imshow(s, hold=False, extent=(yl, yr, xl, xr), origin='lower')
    colorbar(h, orientation='horizontal')
    title('bz, tstep = ' + repr(tstep))
    savefig('bz_' + repr(tstep).zfill(6) + '.png')
        
    # p contour
    figure()
    s = p[:,:,sl]
    h = imshow(s, hold=False, extent=(yl, yr, xl, xr), origin='lower')
    colorbar(h, orientation='horizontal')
    title('p, tstep = ' + repr(tstep))
    savefig('p_' + repr(tstep).zfill(6) + '.png')
        
    # dB contour
    figure()
    s = log10(abs(dB[:,:,sl]))
    h = imshow(s, hold=False, extent=(yl, yr, xl, xr), origin='lower')
    colorbar(h, orientation='horizontal')
    title('log(div B), tstep = ' + repr(tstep))
    savefig('divB_' + repr(tstep).zfill(6) + '.png')
        
    # Jcur contour
    figure()
    s = j[slx,:,:]
    h = imshow(s, hold=False, extent=(zl, zr, xl, xr), origin='lower')
    colorbar(h, orientation='horizontal')
    title('Jcurrent, tstep = ' + repr(tstep))
    savefig('Jcur_' + repr(tstep).zfill(6) + '.png')
        
    # ptot contour
    figure()
    s = ptot[:,:,sl]
    h = imshow(s, hold=False, extent=(yl, yr, xl, xr), origin='lower')
    colorbar(h, orientation='horizontal')
    title('total pressure, tstep = ' + repr(tstep))
    savefig('ptot_' + repr(tstep).zfill(6) + '.png')
    

##### end of script #####
