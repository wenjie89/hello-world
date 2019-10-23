%% analytical solution

close all; clear; clc;

% define parameters

eps0 = 8.85e-12; % vacuum permittivity

Lz = 9e-2; % Lx = 3cm
Lr = 3e-2; % Ly = 3cm

z0 = Lz/10;  % z coordinate of the center of the source sphere
r0 = Lr/20;  % radius of the source sphere

rou0 = 8.85e-12;  % charge density in the sphere

nr = 300;  % number of points 
%nz = nr;


% discretization domain

delr = Lr/200;
delz = delr;

z = 0:delz:Lz;
r = 0:delr:Lr;



[R, Z] = meshgrid(r,z);

%% Analytical solution

idx1 = find(z>0 & z<z0-r0);
idx1 = [find(z==0), idx1];

idx2 = find(z>z0-r0 & z<z0);
idx2 = [find(z==z0-r0), idx2];

idx3 = find(z>z0 & z<z0+r0);
idx3 = [find(z==z0), idx3];

idx4 = find(z>z0+r0 | z==z0+r0);

E_ana = zeros(1,length(z));

E_ana(idx1) = -rou0*r0^3/(3*eps0)./abs(z(idx1)-z0).^2;
E_ana(idx4) =  rou0*r0^3/(3*eps0)./abs(z(idx4)-z0).^2;
E_ana(idx2) = -rou0/(3*eps0)*abs(z(idx2)-z0);
E_ana(idx3) =  rou0/(3*eps0)*abs(z(idx3)-z0);

Eimage = -rou0*r0^3/(3*eps0)./abs(z+z0).^2;

E_ana = E_ana+Eimage;


v_ana = zeros(1,length(z));

v_ana(idx1) = -rou0*r0^3/(3*eps0)*(1./(z(idx1)-z0)+1./(z(idx1)+z0));

v1 = -rou0*r0^3/(3*eps0)*(-1/r0+1/(2*z0-r0));
v_ana(idx2) = v1+rou0/(6*eps0)*r0^2+rou0*r0^3/(3*eps0)/(2*z0-r0)-( rou0/(6*eps0)*(z(idx2)-z0).^2 + rou0*r0^3/(3*eps0)./(z(idx2)+z0) );

v2 = rou0/(6*eps0)*r0^2+rou0*r0^3/(3*eps0)/(2*z0-r0)-rou0*r0^3/(3*eps0)./(2*z0);

v_ana(idx3) = v1+v2 + rou0*r0^3/(3*eps0)/(2*z0) -(rou0/(6*eps0)*(z(idx3)-z0).^2+rou0*r0^3/(3*eps0)./(z(idx3)+z0));

v3 = rou0*r0^3/(3*eps0)/(2*z0)-(rou0/(6*eps0)*r0^2+rou0*r0^3/(3*eps0)/(2*z0+r0));

v_ana(idx4) = v1+v2+v3+rou0*r0^3/(3*eps0)*(-1/r0+1/(2*z0+r0)+1./(z(idx4)-z0)-1./(z(idx4)+z0));


v_ana2 = cumsum(-E_ana*(z(2)-z(1)));





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
xlabel('r (m)')
ylabel('z (m)')
title('electron density distribution')
set(gca, 'Fontsize', 15)

v = zeros(size(R,1),size(R,2));
dnorms = sum(sum(S.^2));
epsilon = 1e-10;
w = 1.9;


% match the boundary condition------------------------------------------------------------

% v(r=0)
v(:,1) = v_ana.';


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
plot(z,v_ana,'r--')
plot(z,v_ana2,'g--')
xlabel('z (m)')
ylabel('electric potential (V)')
legend('Numerical solution', 'Analytical solution','Analytical solution2')
title('electric potential as a function of z at r = 0')
set(gca, 'Fontsize', 15)


E_ana(1:find(z==z0))=-E_ana(1:find(z==z0));


figure(6)
clf
plot(z,E(:,find(r==0)))
hold on 
plot(z,E_ana,'r--')
plot([z0 z0],[0 6]*1e-4,'k-.')
xlim([0 0.05])
xlabel('z (m)')
ylabel('electric field (V/m)')
legend('Numerical solution', 'Analytical solution','z=z0')
title('electric field as a function of z at r = 0')
set(gca, 'Fontsize', 15)





