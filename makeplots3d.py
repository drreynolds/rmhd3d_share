# matplotlib-based plotting script for RMHD3d tests
# Daniel R. Reynolds, reynolds@smu.edu

# imports
from pylab import *
import shlex
import mhd2d_tools as mhd2d

# load the relevant files, and parameters
times, energy, mass = mhd2d.load_cons()
kenergy, menergy, ienergy = mhd2d.load_energies()
xl, xr, yl, yr, zl, zr, ndump = mhd2d.load_mhd()
mu, eta, kappa, gamma = load_props()
hprec, vprec, combo, nwaves, order, ADIcorr = load_prec()
nx, ny, nz, xprocs, yprocs, zprocs, xbc, ybc, zbc = load_mesh()
Newton, Krylov, Fnorm, dB, SolTime = load_RunHistory()

# plot conservation of energy
figure()
plot(times,energy,'b-')
xlabel('time')
ylabel('total energy')
title('Conservation of total energy')
grid()
savefig('energy_cons.png')

# plot conservation of mass
figure()
plot(times,mass,'b-')
xlabel('time')
ylabel('mass')
title('Conservation of mass')
grid()
savefig('mass_cons.png')

# plot energy partitioning
figure()
plot(times,kenergy,'b-',times,menergy,'r--',times,ienergy,'k:')
xlabel('time')
ylabel('energies')
title('Energy partitioning')
grid()
savefig('energies.png')

# plot Newton solver history
figure()
plot(Newton,'b-')
xlabel('Time step')
ylabel('Iterations')
title('Newton Iteration History')
grid()
savefig('newton_hist.png')

# plot Krylov solver history
figure()
plot(Krylov,'b-')
xlabel('Time step')
ylabel('Iterations')
title('Krylov Iteration History')
grid()
savefig('krylov_hist.png')

# plot Fnorm history
figure()
plot(Fnorm,'b-')
xlabel('Time step')
ylabel('||f||')
title('Nonlinear Residual History')
grid()
savefig('fnorm_hist.png')

# plot divB history
figure()
plot(dB,'b-')
xlabel('Time step')
ylabel('max(div B)')
title('Div(B) History')
grid()
savefig('divB_hist.png')

# set the total number of snapshots
Nt = times.size

# only print out ~10 snapshots of each field
if (Nt/ndump > 10):
    ndump *= 2
if (Nt/ndump > 10):
    ndump *= 2
if (Nt/ndump > 10):
    ndump *= 2

# slice through the center of the domain (in z-direction)
sl = int(nz/2)

# loop over snapshots, loading values and plotting
for tstep in range(0,Nt+1,ndump):
    
    # load solution data
    x,y,z,rho,u,v,w,bx,by,bz,p,dB,j,te = mhd2d.load_vals(tstep)

    # rho contour
    figure()
    s = rho[sl][:][:]
    h = imshow(s, hold=False, extent=(xl, xr, yl, yr), origin='lower')
    colorbar(h, orientation='horizontal')
    title('rho, t = ' + repr(round(times[tstep],4)))
    savefig('rho_' + repr(tstep).zfill(6) + '.png')
        
    # u contour
    figure()
    s = u[sl][:][:]
    h = imshow(s, hold=False, extent=(xl, xr, yl, yr), origin='lower')
    colorbar(h, orientation='horizontal')
    title('u, t = ' + repr(round(times[tstep],4)))
    savefig('u_' + repr(tstep).zfill(6) + '.png')
        
    # v contour
    figure()
    s = v[sl][:][:]
    h = imshow(s, hold=False, extent=(xl, xr, yl, yr), origin='lower')
    colorbar(h, orientation='horizontal')
    title('v, t = ' + repr(round(times[tstep],4)))
    savefig('v_' + repr(tstep).zfill(6) + '.png')
        
    # w contour
    figure()
    s = w[sl][:][:]
    h = imshow(s, hold=False, extent=(xl, xr, yl, yr), origin='lower')
    colorbar(h, orientation='horizontal')
    title('w, t = ' + repr(round(times[tstep],4)))
    savefig('w_' + repr(tstep).zfill(6) + '.png')
        
    # bx contour
    figure()
    s = bx[sl][:][:]
    h = imshow(s, hold=False, extent=(xl, xr, yl, yr), origin='lower')
    colorbar(h, orientation='horizontal')
    title('bx, t = ' + repr(round(times[tstep],4)))
    savefig('bx_' + repr(tstep).zfill(6) + '.png')
        
    # by contour
    figure()
    s = by[sl][:][:]
    h = imshow(s, hold=False, extent=(xl, xr, yl, yr), origin='lower')
    colorbar(h, orientation='horizontal')
    title('by, t = ' + repr(round(times[tstep],4)))
    savefig('by_' + repr(tstep).zfill(6) + '.png')
        
    # bz contour
    figure()
    s = bz[sl][:][:]
    h = imshow(s, hold=False, extent=(xl, xr, yl, yr), origin='lower')
    colorbar(h, orientation='horizontal')
    title('bz, t = ' + repr(round(times[tstep],4)))
    savefig('bz_' + repr(tstep).zfill(6) + '.png')
        
    # p contour
    figure()
    s = p[sl][:][:]
    h = imshow(s, hold=False, extent=(xl, xr, yl, yr), origin='lower')
    colorbar(h, orientation='horizontal')
    title('p, t = ' + repr(round(times[tstep],4)))
    savefig('p_' + repr(tstep).zfill(6) + '.png')
        
    # dB contour
    figure()
    s = log10(abs(dB[sl][:][:]))
    h = imshow(s, hold=False, extent=(xl, xr, yl, yr), origin='lower')
    colorbar(h, orientation='horizontal')
    title('log(div B), t = ' + repr(round(times[tstep],4)))
    savefig('divB_' + repr(tstep).zfill(6) + '.png')
        
    # Jcur contour
    figure()
    s = j[sl][:][:]
    h = imshow(s, hold=False, extent=(xl, xr, yl, yr), origin='lower')
    colorbar(h, orientation='horizontal')
    title('Jcurrent, t = ' + repr(round(times[tstep],4)))
    savefig('Jcur_' + repr(tstep).zfill(6) + '.png')
        
    # te contour
    figure()
    s = te[sl][:][:]
    h = imshow(s, hold=False, extent=(xl, xr, yl, yr), origin='lower')
    colorbar(h, orientation='horizontal')
    title('total energy, t = ' + repr(round(times[tstep],4)))
    savefig('te_' + repr(tstep).zfill(6) + '.png')
    

##### end of script #####
