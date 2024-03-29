%% analytical solution

close all; clear; clc;

% define parameters

eps0 = 8.85e-12; % vacuum permittivity

Lz = 9e-2; % Lx = 3cm
Lr = 3e-2; % Ly = 3cm

z0 = Lz/2;  % z coordinate of the center of the source sphere
r0 = Lr/6;  % radius of the source sphere

rou0 = 8.85e-12;  % charge density in the sphere


nr = 100;  % number of points 
nz = nr;


% discretization domain

delr = Lr/200;
delz = delr;

z = 0:delz:Lz;
r = delr/2:delr/2:Lr;




[Z, R] = meshgrid(z,r);

% Compute the electric field and electric potential

% r  = sqrt((Z-x0).^2+(R-y0).^2);


% E1 = ro_0 / (2*eps0) * r;       % r < R
% E2 = ro_0*R^2/(2*eps0) ./ r;    % r > R
% 
% Phi1 = -ro_0/(4*eps0) * r.*2;   % r < R
% Phi2 = -ro_0/(4*eps0) * R^2 * (1-2*log(R./r)); % r > R 


% analytical solution 

% E1 = rou0 / (2*eps0) * sqrt((Z-x0).^2+(R-y0).^2);       % r < R
% E2 = rou0*R^2/(2*eps0) ./ sqrt((Z-x0).^2+(R-y0).^2);    % r > R
% 
% V1 = -rou0/(4*eps0) * ((Z-x0).^2+(R-y0).^2);   % r < R
% V2 = -rou0/(4*eps0) * R^2 * (1-2*log(R./sqrt((Z-x0).^2+(R-y0).^2))); % r > R 
% 
% [idxX, idxY] = find( sqrt( (Z-x0).^2 + (R-y0).^2 ) > R );
% 
% 
% E = E1;
% V = V1;
% 
% for i = 1:length(idxX) 
%     E(idxX(i), idxY(i)) = E2(idxX(i), idxY(i));
%     V(idxX(i), idxY(i)) = V2(idxX(i), idxY(i));
% end
% 
% 
% Ex = E.* (Z-x0)./r;   % Ex = Er * (x-x0)/r
% Ey = E.* (R-y0)./r;   % Ey = Er * (y-y0)/r

% xidx = find(r == y0);
% yidx = find(z == x0);

%% numerical solution

% define the source sphere 
S = (delr*delz)^2/eps0/(delr^2+delz^2)/2 * rou0 * ones(size(Z,1),size(Z,2));

sourceRegion = (Z-z0).^2 + R.^2 - r0^2; 

[idxX, idxY] = find(sourceRegion > 0);

for i = 1:length(idxX)
    S(idxX(i),idxY(i)) = 0;
end

% plot the source region
figure(1)
clf
imagesc(r,z,S.')
colorbar
xlabel('r')
ylabel('z')
title('electron density distribution')
set(gca, 'Fontsize', 15)
% xlim([0 r0])
% axis equal

% define coefficients for the SOR method
eleCoe = -(delr*delz)^2/2/(delz^2+delr^2);


a = zeros(1,length(r));

for j = 1:length(r)
    a(j) = eleCoe*(-1/(2*r(j)*delr) + 1/delr^2);
    b(j) = eleCoe*( 1/(2*r(j)*delr) + 1/delr^2);
end

c = -delr^2/(delz^2+delr^2)/2;
d = c;
v = zeros(size(R,1),size(R,2));


dnorms = sum(sum(S.^2));
epsilon = 1e-10;
w = 1.9;


[v,it,dnorm] = sorEE430S19(a,b,c,d,v,S,w,epsilon,nz,nr,dnorms);


dnorm=0;
for iz = 2:nz
    for jr = 2:nr

        residual = a(jr)*v(iz,jr-1)+b(jr)*v(iz,jr+1)+...
                c*v(iz-1,jr)+d*v(iz+1,jr)+v(iz,jr)-S(iz,jr);
        dnorm=dnorm+residual^2;
    end
end









