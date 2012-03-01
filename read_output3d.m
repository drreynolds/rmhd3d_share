function [rho,u,v,w,Bx,By,Bz,p,divB,J,e,x,y,z] = read_output3d(fname,Nx,Ny,Nz)
%
% [rho,u,v,w,Bx,By,Bz,p,divB,J,e,x,y,z] = read_output3d(fname,Nx,Ny,Nz)
%
% Extracts 3D fields from gnuplot output file.
%
% Inputs:  
%     fname - file name for input
%        Nx - dimension of fields in x-direction
%        Ny - dimension of fields in y-direction
%        Nz - dimension of fields in z-direction
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
%         e - total energy field (3D)
%         x - x-coordinates (1D)
%         y - y-coordinates (1D)
%         z - z-coordinates (1D)
%
% Daniel Reynolds
% Mathematics @ SMU
%

% open file for input
fin = fopen(fname,'r');
if fin < 0
   error([ 'Could not open ',fname,' for input']);
end

% check whether there is a label line (marked by #)
nhead = 0;
buffer = fgetl(fin);
[text,buffer] = strtok(buffer);
if (text == '#') 
   nhead = 1;
end
frewind(fin);

% discard header text in preparation of reading data
for i=1:nhead,  buffer=fgetl(fin);  end

% initialize output arrays
rho = zeros(Nx,Ny,Nz);
u = zeros(Nx,Ny,Nz);
v = zeros(Nx,Ny,Nz);
w = zeros(Nx,Ny,Nz);
Bx = zeros(Nx,Ny,Nz);
By = zeros(Nx,Ny,Nz);
Bz = zeros(Nx,Ny,Nz);
p = zeros(Nx,Ny,Nz);
divB = zeros(Nx,Ny,Nz);
J = zeros(Nx,Ny,Nz);
e = zeros(Nx,Ny,Nz);
x = zeros(Nx,1);
y = zeros(Ny,1);
z = zeros(Nz,1);

% read in data
for k=1:Nz
   for j=1:Ny
      for i=1:Nx
	 
	 % read in line of text
	 buffer = fgetl(fin);
	 
	 % extract fields into appropriate output arrays
	 [x(i), y(j), z(k), rho(i,j,k), u(i,j,k), v(i,j,k), w(i,j,k), Bx(i,j,k), ...
	  By(i,j,k), Bz(i,j,k), p(i,j,k), divB(i,j,k), J(i,j,k), e(i,j,k)] ...
          = strread(buffer,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f');
	 
	 % note: gnuplot output designed for 2D, so no z coordinate is saved
	 z(k) = k;
      end
      
      % skip a line
      buffer = fgetl(fin);
   
   end
end


% close file
fclose(fin);
