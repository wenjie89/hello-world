%% analytical solution

close all; clear; clc;

% define parameters

eps0 = 8.85e-12; % vacuum permittivity

Lz = 9e-2; % Lx = 3cm
Lr = 3e-2; % Ly = 3cm

z0 = Lz/2;  % z coordinate of the center of the source sphere
r0 = Lr/6;  % radius of the source sphere
hi
rou0 = 8.85e-12;  % charge density in the sphere

nr = 100;  % number of points 
nz = nr;


% discretization domain

delr = Lr/200;
delz = delr;

z = 0:delz:Lz;
r = 0:delr:Lr;



[R, Z] = meshgrid(r,z);

%% Analytical solution
idxz1 = find(abs(z-z0)<r0 | abs(z-z0)==r0);
idxz2 = find(abs(z-z0)>r0);

V_ana_z = zeros(1,length(z));
V_ana_z(idxz1) = -rou0./(6*eps0)*(z(idxz1)-z0).^2;
V_ana_z(idxz2) = -rou0*r0^2/(6*eps0)*(3-2*r0./abs(z(idxz2)-z0));

E_ana_z = zeros(1,length(z));
E_ana_z(idxz1) = rou0/(3*eps0).*abs(z(idxz1)-z0);
E_ana_z(idxz2) = rou0*r0^3./(3*eps0*(z(idxz2)-z0).^2);



idxr1 = find(r<r0 | r==r0);
idxr2 = find(r>r0);

V_ana_r = zeros(1,length(r));
V_ana_r(idxr1) = -rou0./(6*eps0)*r(idxr1).^2;
V_ana_r(idxr2) = -rou0*r0^2/(6*eps0)*(3-2*r0./r(idxr2));

E_ana_r  = zeros(1,length(r));
E_ana_r(idxr1) = rou0/(3*eps0).*r(idxr1);
E_ana_r(idxr2) = rou0*r0^3./(3*eps0*r(idxr2).^2);



%% numerical solution


% define the source sphere  (it's good to test the cylinder first)------------------------
S = rou0 * ones(size(Z,1),size(Z,2));

sourceRegion = (Z-z0).^2 + R.^2 - r0^2; 

[idxX, idxY] = find(sourceRegion > 0);

for i = 1:length(idxX)
    S(idxX(i),idxY(i)) = 0;
end
%S = S';

% plot the source region
figure(1)
clf
imagesc(r,z,S)
colorbar
xlabel('r')
ylabel('z')
title('electron density distribution')
set(gca, 'Fontsize', 15)


v = zeros(size(R,1),size(R,2));
dnorms = sum(sum(S.^2));
epsilon = 1e-10;
w = 1.9;


% match the boundary condition------------------------------------------------------------

% v(r=0)
v(idxz1,1) = -rou0/(6*eps0)*(z(idxz1)-z0).^2;
v(idxz2,1) = -rou0*r0^2/(6*eps0)*(3-2*r0./abs(z(idxz2)-z0));

% v(r=Lr)  assume Lr>r0
v(:,size(v,2)) = -rou0 * r0^2/(6*eps0) * (3-2*r0./ sqrt( (z-z0).^2+Lr^2 ) ); 

% v(z=0) assume z0>r0
v(1,:) = -rou0 * r0^2/(6*eps0) * (3-2*r0./ sqrt( z0^2+r.^2 ) );

% v(z=Lz) assume Lz>z0+r0
v(size(v,1),:) = -rou0 * r0^2/(6*eps0) * (3-2*r0./ sqrt( (Lz-z0)^2+r.^2 ) );


% SOR-------------------------------------------------------------------------------------
[v,it] = mySOR_cylinder(v,r,z,S,w,epsilon);
[Er, Ez] = gradient(v);
Er = -Er/delr;
Ez = -Ez/delz;
E = sqrt(Er.^2+Ez.^2);


%% figures

figure(2)
clf
imagesc(r,z,v)
colorbar
xlabel('r (m)')
ylabel('z (m)')
title('electric potential (V)')
set(gca, 'Fontsize', 15)

figure(3)
clf
imagesc(r,z,E)
colorbar
xlabel('r (m)')
ylabel('z (m)')
title('electric field (V/m)')
set(gca, 'Fontsize', 15)


figure(4)
clf
plot(z,v(:,find(r==0)))
hold on 
plot(z,V_ana_z,'r--')
xlabel('z (m)')
ylabel('electric potential (V)')
legend('Numerical solution', 'Analytical solution')
title('electric potential as a function of z at r = 0')
set(gca, 'Fontsize', 15)

figure(5)
clf
plot(r,v(find(z==z0),:))
hold on 
plot(r,V_ana_r,'r--')
xlabel('r (m)')
ylabel('electric potential (V)')
legend('Numerical solution', 'Analytical solution')
title('electric potential as a function of r at z = z_0')
set(gca, 'Fontsize', 15)


figure(6)
clf
plot(z,E(:,find(r==0)))
hold on 
plot(z,E_ana_z,'r--')
xlabel('z (m)')
ylabel('electric field (V/m)')
legend('Numerical solution', 'Analytical solution')
title('electric field as a function of z at r = 0')
set(gca, 'Fontsize', 15)

figure(7)
clf
plot(r,E(find(z==z0),:))
hold on 
plot(r,E_ana_r,'r--')
xlabel('r (m)')
ylabel('electric field (V/m)')
legend('Numerical solution', 'Analytical solution')
title('electric field as a function of r at z = z_0')
set(gca, 'Fontsize', 15)



