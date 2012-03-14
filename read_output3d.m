function [rho,u,v,w,Bx,By,Bz,p,divB,J,ptot,x,y,z] = read_output3d(idump)
%
% [rho,u,v,w,Bx,By,Bz,p,divB,J,ptot,x,y,z] = read_output3d(idump)
%
% Extracts 3D fields from gnuplot output file.
%
% Inputs:  
%     idump - output number to load (e.g. data000002 would use "2")
%
% Outputs:  
%       rho - density field (3D)
%         u - x-velocity field (3D)
%         v - y-velocity field (3D)
%         w - z-velocity field (3D)
%        Bx - x-magnetic field (3D)
%        By - y-magnetic field (3D)
%        Bz - z-magnetic field (3D)
%         p - pressure field (3D)
%      divB - div(B) field [for constraint check] (3D)
%         J - toroidal current (3D)
%      ptot - total pressure field (3D)
%         x - x-coordinates (1D)
%         y - y-coordinates (1D)
%         z - z-coordinates (1D)
%
% Daniel Reynolds
% Mathematics @ SMU
%

% convert idump to determine the data directory
if (idump<10)
   sdump = [ 'data00000' , sprintf('%i',idump) , '/' ];
elseif (idump<100)
   sdump = [ 'data0000' ,  sprintf('%i',idump) , '/' ];
elseif (idump<1000)
   sdump = [ 'data000' ,   sprintf('%i',idump) , '/' ];
elseif (idump<10000)
   sdump = [ 'data00' ,    sprintf('%i',idump) , '/' ];
elseif (idump<100000)
   sdump = [ 'data0' ,     sprintf('%i',idump) , '/' ];
else
   sdump = [ 'data' ,      sprintf('%i',idump) , '/' ];
end

% load input files
%   get processor layout
proclayout = load('proc_layout.txt');
%   extract xprocs, yprocs, zprocs
xprocs = proclayout(1,2);
yprocs = proclayout(1,3);
zprocs = proclayout(1,4);
%   get x grid
x = load('x.grid.001');
for i=2:xprocs
   if (i<10)
      s = [ 'x.grid.00' , sprintf('%i',i) ];
   elseif (i<100)
      s = [ 'x.grid.0' ,  sprintf('%i',i) ];
   else
      s = [ 'x.grid.' ,   sprintf('%i',i) ];
   end
   x = [x; load(s)];
end
%   get y grid
y = load('y.grid.001');
for i=2:yprocs
   if (i<10)
      s = [ 'y.grid.00' , sprintf('%i',i) ];
   elseif (i<100)
      s = [ 'y.grid.0' ,  sprintf('%i',i) ];
   else
      s = [ 'y.grid.' ,   sprintf('%i',i) ];
   end
   y = [y; load(s)];
end
%   get z grid
z = load('z.grid.001');
for i=2:zprocs
   if (i<10)
      s = [ 'z.grid.00' , sprintf('%i',i) ];
   elseif (i<100)
      s = [ 'z.grid.0' ,  sprintf('%i',i) ];
   else
      s = [ 'z.grid.' ,   sprintf('%i',i) ];
   end
   z = [z; load(s)];
end

% get overall dimensions
Nx = length(x);
Ny = length(y);
Nz = length(z);

% set Nxlocal, Nylocal, Nzlocal as local problem sizes (general)
Nxlocal = Nx/xprocs;
Nylocal = Ny/yprocs;
Nzlocal = Nz/zprocs;

% initialize output arrays
rho  = zeros(Nx,Ny,Nz);
u    = zeros(Nx,Ny,Nz);
v    = zeros(Nx,Ny,Nz);
w    = zeros(Nx,Ny,Nz);
Bx   = zeros(Nx,Ny,Nz);
By   = zeros(Nx,Ny,Nz);
Bz   = zeros(Nx,Ny,Nz);
p    = zeros(Nx,Ny,Nz);
divB = zeros(Nx,Ny,Nz);
J    = zeros(Nx,Ny,Nz);
ptot = zeros(Nx,Ny,Nz);

% number of fields to extract
nvar = 11;

% iterate over output files
for i=0:xprocs*yprocs*zprocs-1
   % choose output file
   if (i<10)
      fname = [ sdump, 'output00000' , sprintf('%i',i) , '.txt' ];
   elseif (i<100)
      fname = [ sdump, 'output0000' ,  sprintf('%i',i) , '.txt' ];
   elseif (i<1000)
      fname = [ sdump, 'output000' ,   sprintf('%i',i) , '.txt' ];
   elseif (i<10000)
      fname = [ sdump, 'output00' ,    sprintf('%i',i) , '.txt' ];
   elseif (i<100000)
      fname = [ sdump, 'output0' ,     sprintf('%i',i) , '.txt' ];
   else
      fname = [ sdump, 'output' ,      sprintf('%i',i) , '.txt' ];
   end

   % load output file
   data = load(fname);

   % determine processor location in overall grid
   iprocx = proclayout(i+2,2);
   iprocy = proclayout(i+2,3);
   iprocz = proclayout(i+2,4);

   % determine local processor extents in global grid
   is = (iprocx-1)*Nxlocal + 1;
   ie = min(is+Nxlocal-1, Nx);
   js = (iprocy-1)*Nylocal + 1;
   je = min(js+Nylocal-1, Ny);
   ks = (iprocz-1)*Nzlocal + 1;
   ke = min(ks+Nzlocal-1, Nz);

   % determine local problem sizes (and check for consistency)
   myNx = ie-is+1;
   myNy = je-js+1;
   myNz = ke-ks+1;
   N = length(data);
   if (N ~= myNx*myNy*myNz*nvar)
      disp(sprintf('read_output: incompatible file dims p%i: (%i,%i,%i,%i) vs %i',...
          i,myNx,myNy,myNz,nvar,N));
   end

   % reshape local data into 4D array
   data = reshape(data,myNx,myNy,myNz,nvar);

   % insert results into output arrays
   rho(is:ie,js:je,ks:ke)  = data(:,:,:,1);
   u(is:ie,js:je,ks:ke)    = data(:,:,:,2);
   v(is:ie,js:je,ks:ke)    = data(:,:,:,3);
   w(is:ie,js:je,ks:ke)    = data(:,:,:,4);
   bx(is:ie,js:je,ks:ke)   = data(:,:,:,5);
   by(is:ie,js:je,ks:ke)   = data(:,:,:,6);
   bz(is:ie,js:je,ks:ke)   = data(:,:,:,7);
   p(is:ie,js:je,ks:ke)    = data(:,:,:,8);
   divB(is:ie,js:je,ks:ke) = data(:,:,:,9);
   J(is:ie,js:je,ks:ke)    = data(:,:,:,10);
   ptot(is:ie,js:je,ks:ke) = data(:,:,:,11);

end
