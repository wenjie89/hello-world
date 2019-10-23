%%
close all; clear; clc;

eps0 = 8.85e-12; % vacuum permitivity
rou0 = eps0;


z = 0:0.001:100e-2;   % 0--100 cm

z0 = 3e-2; % m
R = 1.5e-2;

idx1 = find(z>0 & z<z0-R);
idx1 = [find(z==0), idx1];

idx2 = find(z>z0-R & z<z0);
idx2 = [find(z==z0-R), idx2];

idx3 = find(z>z0 & z<z0+R);
idx3 = [find(z==z0), idx3];

idx4 = find(z>z0+R | z==z0+R);


E = zeros(1,length(z));

E(idx1) = -rou0*R^3/(3*eps0)./abs(z(idx1)-z0).^2;
E(idx4) =  rou0*R^3/(3*eps0)./abs(z(idx4)-z0).^2;
E(idx2) = -rou0/(3*eps0)*abs(z(idx2)-z0);
E(idx3) =  rou0/(3*eps0)*abs(z(idx3)-z0);

Eimage = -rou0*R^3/(3*eps0)./abs(z+z0).^2;

E = E+Eimage;

phi = cumsum(-E*(z(2)-z(1)));





%%
figure(1)
clf
plot(z,E)
hold on 
plot([z0-R, z0-R],[-6 6]*1e-3)
plot([z0+R, z0+R],[-6 6]*1e-3)

figure(2)
clf
plot(z,phi)
hold on
plot([z0-R, z0-R], [0 1]*1e-4)
plot([z0+R, z0+R] ,[0 1]*1e-4)



















