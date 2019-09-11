close all; clear; clc;

%% define parameters

eps0 = 8.85e-12; % vacuum permittivity
Lx = 3e-2; % Lx = 3cm
Ly = 3e-2; % Ly = 3cm

R = Lx/6;  % radius of the cylinder
ro_0 = 8.85e-12;  % charge density in the cylinder


nx = 20;  % number of points 
ny = 20;


% discretization domain

x = linspace(0,Lx,nx+1);
y = linspace(0,Ly,ny+1);

[X, Y] = meshgrid(x,y);

%% Compute the electric field and electric potential


E1 = ro_0/(2*eps0) * sqrt((X-Lx/2).^2+(Y-Ly/2).^2);       % r < R
E2 = ro_0*R^2/(2*eps0) ./ sqrt((X-Lx/2).^2+(Y-Ly/2).^2);  % r > R

Phi1 = -ro_0/(4*eps0) * ((X-Lx/2).^2+(Y-Ly/2).^2);  % r < R
Phi2 = -ro_0/(4*eps0) * ((X-Lx/2).^2+(Y-Ly/2).^2).*...
      (1+2*log( sqrt((X-Lx/2).^2+(Y-Ly/2).^2)/R )); % r > R 

[idxX, idxY] = find( sqrt( (X-Lx/2).^2 + (Y-Ly/2).^2 ) > R );


E = E1;
Phi = Phi1;

for i = 1:length(idxX) 
    E(idxX(i), idxY(i)) = E2(idxX(i), idxY(i));
    Phi(idxX(i), idxY(i)) = Phi2(idxX(i), idxY(i));
end


Ex = E.* (X-Lx/2)/R;
Ey = E.* (Y-Ly/2)/R;




%%

figure(1)
clf
contour(X,Y,Phi)
hold on
quiver(X,Y,Ex,Ey)
axis equal

figure(2)
clf
plot(x,E(find(y==Ly/2),:))

figure(3)
clf
plot(x,Phi(find(y==Ly/2),:))





























