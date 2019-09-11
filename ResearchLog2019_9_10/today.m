%% define basic parameters

close all; clear; clc;



nr = 151;                         % no. of points in r direction
nz = 1681;                        % no. of points in z direction
zmin = 0; zmax = 1.4e-2;          % computational domain in z direction [m]
rmin = -0.125e-2; rmax = 0.125e-2;        % computational domain in r direction [m]
z0 = 0.7e-2;                      % middle of the compurational domain [m]
n0 = 1e20;                        % number density of electrons


sigma_r = 0.02e-2;  
sigma_z = 0.02e-2;


vz = 1;                          % velocity in z direction
vr = 1;                          % velocity in r direction
tmin = 0; 
tmax = 3e-4;
dt = 5e-6;                       % time resolution





%% discretize computational domain


rVec = linspace(rmin,rmax,nr);
zVec = linspace(zmin,zmax,nz);  
[r, z] = meshgrid(rVec,zVec);
nTimeStep = floor((tmax-tmin)/dt);


%% build the initial function (the number density of electrons) to be transported

ne = n0 * exp( -(r/sigma_r).^2 - ( (z-z0) / sigma_z ).^2 );

figure(1)
imagesc(rVec,zVec,ne)
colorbar
xlabel('(m)')
ylabel('(m)')
title('density distribution of the electrons')

figure(2)
plot(z(:,find(rVec==0)),ne(:,find(rVec==0)))
xlabel('(m)')
ylabel('electron density')

figure(3)
plot(r(find(zVec==(zmin+zmax)/2),:), ne(find(zVec==(zmin+zmax)/2),:))
xlabel('(m)')
ylabel('electron density')



%% main algorithm

% for n = 1:nTimeStep
%     
%     for i = 1:length(rVec)
%         
%         for j = 1:length(zVec)


% step 1 -- low order scheme






% step 2 -- high order scheme


% step 3 -- antidiffusion


% step 4 -- compute updated low order solution


% step 5 -- limiter 


% step 6 -- limited antidiffusive fluxes


















