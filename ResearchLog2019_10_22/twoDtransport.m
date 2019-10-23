close all; clear; clc;

delete *.avi 

% define parameters

 

eps0 = 8.85e-12; % vacuum permittivity

q = 1.602e-19;   % [C] elementary charge

 

Lz = 9e-2; % Lx = 3cm

Lr = 3e-2; % Ly = 3cm

 

z0 = Lz/2;  % z coordinate of the center of the source sphere

r0 = Lr/6;  % radius of the source sphere

 

rou0 = 8.85e-12;  % [C/m^3]  charge density in the sphere

 

nr = 100;  % number of points 

nz = nr;



% domain  discretization 

 

delr = Lr/200;

delz = delr;

 

z = 0:delz:Lz;

r = 0:delr:Lr;

 


[R, Z] = meshgrid(r,z);

 


%% numerical solution

 

 

% define the source sphere  (it's good to test the cylinder first)------------------------

S = rou0 * ones(size(Z,1),size(Z,2));  %  [C/m^3]

 

sourceRegion = (Z-z0).^2 + R.^2 - r0^2; 

 

[idxX, idxY] = find(sourceRegion > 0);

 

for i = 1:length(idxX)

    S(idxX(i),idxY(i)) = 0;

end

 

ne = S/q;   % [m^-3] density distribution of positive charges

 

% plot the source region

figure(1)

clf

imagesc(r,z,ne)

colorbar

xlabel('r')

ylabel('z')

title('electron density distribution')

set(gca, 'Fontsize', 15)

 

 

%% the main algorithm

 

% define parameters --------------------------------------------------

 

vz = 1; % [m/s]

vr = 0.0005; 

dt = 0.05; % [x]

% dx = 1e-2; % [m]

t = 0; % current time [s]

T = 10; % maximum time [s]

 

 

 

% transport the pulse

 

ntimestep = T/dt;

 

 

z = 0:delz:Lz; % z0.5, z1.5, z2.5...

r = 0:delr:Lr; % r1, r2, r3...

 


Sz = zeros(1, size(ne,2)); 

Sr = zeros(1, size(ne,2)); 

for j = 2:size(ne,2)

    Sz(j) = pi*((r(j)+delr/2)^2-(r(j)-delr/2)^2);  % Sz2.5, Sz3.5...

    Sr(j) = 2*pi*(r(j)+delr/2)*delz;               % Sr2.5, Sr3.5...

end

Sz(1) = pi*(delr/2)^2;       % Sz1.5

Sr(1) = 2*pi*(delr/2)*delz;  % Sr1.5

 

Sz = repmat(Sz,size(ne,1),1);

Sr = repmat(Sr,size(ne,1),1);

V = Sz*delz;

 

 

%% Set up the movie.
writerObj = VideoWriter('out.avi'); % Name it.
writerObj.FrameRate = 60; % How many frames per second.
open(writerObj);  

% loop through time index

neNew = zeros(size(ne,1), size(ne,2));  % initialize it for every time step

for n = 1:ntimestep 
   

    % numerical solution
 

    for i = 1:size(ne,1)

        neNew(i,1) = ne(i,1) - dt/V(i,1)*( ne(i,1)*vr*Sr(i,1) ); 

    end

    

    for i = 1:size(ne,1)

        for j = 2:size(ne,2)

            neNew(i,j) = ne(i,j) - dt/V(i,j)*(ne(i,j)*vr*Sr(i,j)- ne(i,j-1)*vr*Sr(i,j-1));

        

%         uinp1 = u(i) - v*dt/dx * (u(i) - u(i-1));

%         uinp1Vec(i) = uinp1;

 

        end

    end
    

    ne = neNew;   % update the numerical solution in each loop 

    t = t+dt;     % update the current time index

    

    figure(2)
    


    imagesc(r,z,ne)

    colorbar

    xlabel('r (m)')

    ylabel('z (m)')

    title(sprintf('electron density distribution at t=%.3f s',t))

    set(gca, 'Fontsize', 15)

    %pause(0.1)
    
 %if mod(i,4)==0, % Uncomment to take 1 out of every 4 frames.
        frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
        writeVideo(writerObj, frame);
    %end

    

end
close(writerObj); % Saves the movie.





