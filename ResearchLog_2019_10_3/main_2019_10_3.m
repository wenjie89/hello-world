%% basic settings
close all; clear; clc;

% define parameters 
eps0 = 8.85e-12;  % vacuum permittivity
rou0 = 8.85e-12;  % charge density in the sphere

Lz = 9e-2; % Lz = 9cm
Lr = 3e-2; % Lr = 3cm


nr = 101;  % number of points 
nz = nr;

% domain discretization 

delr = Lr/(nr-1);
delz = Lz/(nz-1);

z = 0:delz/2:Lz;  % z_odd -- interface; z_even -- center
r = 0:delr/2:Lr;  % r_odd -- center; r_even -- interface

% z0 = Lz/2;  % z coordinate of the center of the source sphere
% r0 = Lr/6;  % radius of the source sphere

[R, Z] = meshgrid(r,z);


% define a infinite long cylinder filled with electrons
R0 = Lr/3;  
rou = rou0 * ones(size(R,1), size(R,2));
[IdxI, IdxJ] = find(R>Lr/3);

for i = 1:length(IdxI)
    rou(IdxI(i), IdxJ(i)) = 0;
end

figure(1)
clf
imagesc(r,z,rou)
xlabel('r')
ylabel('z')
title('source distribution')
colorbar






%% define parameters for SOR

u = -2/delr^2-2/delz^2;
a(j) = (-1/2/r(j)/delr + 1/delr^2)/u;
b(j) = ( 1/2/r(j)/delr + 1/delr^2)/u;
c = 1/delz^2/u;
d = c;
s = -rou/eps0/u;



















