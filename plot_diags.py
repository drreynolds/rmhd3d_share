# matplotlib-based plotting script for RMHD3d tests
# Daniel R. Reynolds, reynolds@smu.edu

# imports
from pylab import *
import shlex
import mhd3d_tools as mhd3d

# load the relevant files, and parameters
times, energy, mass = mhd3d.load_cons()
kenergy, menergy, ienergy = mhd3d.load_energies()
Newton, Krylov, Fnorm, dB, SolTime = mhd3d.load_RunHistory()

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


##### end of script #####
