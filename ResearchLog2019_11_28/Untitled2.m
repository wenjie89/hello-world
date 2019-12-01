%% analytical solution

close all; clear; clc;
% tin can

% define parameters

eps0 = 8.85e-12; % vacuum permittivity

Lz = 9e-2; % Lx = 3cm
Lr = 3e-2; % Ly = 3cm

z0 = Lz/2;  % z coordinate of the center of the source sphere
r0 = Lr/6;  % radius of the source sphere
%r0 = 0.9*Lr;  % radius of the source sphere

rou0 = 8.85e-12;  % [C/m^3] charge density in the sphere

nr = 501;  % number of points 



% discretization domain

dr = Lr/(nr-1);
dz = dr;

z = dz/2:dz:Lz;
r = 0:dr:Lr;



[R, Z] = meshgrid(r,z);



%% numerical solution


% define the source sphere  (it's good to test the cylinder first)------------------------
S = rou0 * ones(size(Z,1),size(Z,2));

sourceRegion = (Z-z0).^2 + R.^2 - r0^2; 

[idxX, idxY] = find(sourceRegion > 0);

for i = 1:length(idxX)
    S(idxX(i),idxY(i)) = 0;
end
%S = S';

% volume elements for the whole computational domain ----------------------------
V = zeros(size(Z,1),size(Z,2));  
V(:,1) = pi*(dr/2)^2*dz;
for j=2:size(V,2)
    V(:,j) = pi * ( (r(j)+dr/2 )^2 - ( r(j-1)+dr/2 )^2 )  * dz;
end


% plot the source region
figure(1)
clf
imagesc(r,z,S)
colorbar
xlabel('r (m)')
ylabel('z (m)')
title('electron density distribution (C/m^3)')
set(gca, 'Fontsize', 15)





v = zeros(size(R,1),size(R,2));
dnorms = sum(sum(S.^2));
epsilon = 1e-10;
w = 1.9;


% match the boundary condition---------------------------------------------


% right boundary: v(r=Lr)  assume Lr>r0


for i = 1:size(v,1)
    v(i,size(v,2)) = CalcBoundaryPotential(z(i),Lr,Z,R,S,V);
end

rightBC_ana = rou0*r0^3/(3*eps0)./sqrt( (z-z0).^2+Lr^2 );

figure(2)
clf
plot(z,rightBC_ana)
hold on
plot(z,v(:,size(v,2)),'r--')
legend('analytical', 'numerical')


% bottom boundary: v(z=dz/2) assume z0>r0


for i = 1:size(v,2)
    v(1,i) = CalcBoundaryPotential(0,r(i),Z,R,S,V);
end

bottomBC_ana = rou0 * r0^3/(3*eps0) ./ sqrt( z0^2+r.^2 );

figure(3)
clf
plot(r,bottomBC_ana)
hold on
plot(r,v(1,:),'r--')
legend('analytical', 'numerical')


% top boundary: v(z=Lz) assume Lz>z0+r0


for i = 1:size(v,2)
    v(size(v,1),i) = CalcBoundaryPotential(Lz,r(i),Z,R,S,V);
end

topBC_ana =  rou0 * r0^3/(3*eps0) ./ sqrt( (Lz-z0)^2+r.^2 );

figure(4)
clf
plot(r,topBC_ana)
hold on
plot(r,v(size(v,1),:),'r--')
legend('analytical', 'numerical')



% SOR-------------------------------------------------------------------------------------
[v,it] = mySOR_cylinder(v,r,z,S,w,epsilon);
[Er, Ez] = gradient(v);


Er = -Er/dr;
Ez = -Ez/dz;





