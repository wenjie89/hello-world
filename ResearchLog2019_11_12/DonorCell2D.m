%% 2d donor cell algorighm for the transport equation

close all; clear; clc;


% define parameters

eps0 = 8.85e-12; % vacuum permittivity
q = 1.602e-19;   % [C] elementary charge

Lz = 9e-2; % Lx = 3cm
Lr = 3e-2; % Ly = 3cm

z0 = Lz/2;  % z coordinate of the center of the source sphere
r0 = Lr/6;  % radius of the source sphere

rou0 = 8.85e-12;  % [C/m^3]  charge density in the sphere

nr = 301;  % number of points 
% nz = nr;


% discretization domain

dr = Lr/(nr-1);
dz = dr;

% z = 0:dz:Lz;
% r = 0:dr:Lr;

z = dz/2:dz:Lz; % z1, z2, z3...
r = 0:dr:Lr;    % r1, r2, r3...


[R, Z] = meshgrid(r,z);




%% numerical solution


% define the source sphere  (it's good to test the cylinder first)------------------------
S = rou0 * ones(size(Z,1),size(Z,2));  %  [C/m^3]

sourceRegion = (Z-z0).^2 + R.^2 - r0^2; 

[idxX, idxY] = find(sourceRegion > 0);

for i = 1:length(idxX)
    S(idxX(i),idxY(i)) = 0;
end

N = S/q;   % [m^-3] density distribution of positive charges


% plot the source region
figure(1)
clf
imagesc(r,z,N)
colorbar
xlabel('r (m)')
ylabel('z (m)')
title('charge density distribution')
set(gca, 'Fontsize', 15)



%% the main algorithm

% define parameters --------------------------------------------------

CFLz = 0.9;
%CFLz = 0;
CFLr = 0.5;
%CFLr = 0;

dt = 0.05; % [s]

% vz = CFLz*dz/dt;  % [m/s]
% vr = CFLr*dr/dt; 


t = 0; % current time [s]
T = 15; % maximum time [s]



% transport the pulse

ntimestep = T/dt;


z = dz/2:dz:Lz; % z1, z2, z3...
r = 0:dr:Lr;    % r1, r2, r3...

V  = zeros(size(N,1),size(N,2));                        % initialize V_ij
V(:,1) = pi * ((r(1)+r(2))/2)^2 * dz;                   % V_i1
for j = 2:size(V,2)
    V(:,j) = pi * ((r(j)+dr/2)^2-(r(j)-dr/2)^2) * dz;   % V_ij
end

% figure
% plot(V(1,:))




% assume constant velocity in both directions

vr = CFLr*dr/dt * ones(size(N,1),size(N,2)); 
vz = CFLz*dz/dt * ones(size(N,1),size(N,2)); 
return
Fz = zeros(size(N,1),size(N,2));
Fr = zeros(size(N,1),size(N,2));


% Nnew = zeros(size(N,1),size(N,2));  % initialize Nnew matrix
Nnew = N;  % initialize Nnew matrix


for n = 1:ntimestep
    
    % calculate Fz ---------------------------------------------------------------------
    
    for i = 1:size(N,1)-1
        if vz(i,1)+vz(i+1,1)>0 || vz(i,1)+vz(i+1,1)==0       
            Fz(i,1) = N(i,1) * ( vz(i,1)+vz(i+1,1) )/2 * pi * ((r(1)+r(2))/2)^2;
        else
            Fz(i,1) = N(i+1,1) * ( vz(i,1)+vz(i+1,1) )/2 * pi * ((r(1)+r(2))/2)^2;
        end
    end
    
    for i = 1:size(N,1)-1
        for j = 2:size(N,2)
            if vz(i,j)+vz(i+1,j)>0 || vz(i,j)+vz(i+1,j)==0     
                Fz(i,j) =  N(i,j) * ( vz(i,j)+vz(i+1,j) )/2 * pi * ( (r(j)+dr/2)^2-(r(j)-dr/2)^2 );
            else
                Fz(i,j) =  N(i+1,j) * ( vz(i,j)+vz(i+1,j) )/2 * pi * ( (r(j)+dr/2)^2-(r(j)-dr/2)^2 );
            end
        end
    end
    
    % calculate Fr ---------------------------------------------------------------------
    for i = 1:size(N,1)-1
        if vr(i,1)+vr(i,2)>0 || vr(i,1)+vr(i,2)==0
            Fr(i,1) = N(i,1) * ( vr(i,1)+vr(i,2) ) * pi * dz * (r(1)+r(2))/2;
        else
            Fr(i,1) = N(i,2) * ( vr(i,1)+vr(i,2) ) * pi * dz * (r(1)+r(2))/2;
        end
            
    end
    
    for i = 1:size(N,1)-1
        for j = 2:size(N,2)-1
            if vr(i,j)+vr(i,j+1)>0 || vr(i,j)+vr(i,j+1)==0
                Fr(i,j) = N(i,j) * ( vr(i,j)+vr(i,j+1) ) * pi * dz * (r(j)+r(j+1))/2;
            else
                Fr(i,j) = N(i,j+1) * ( vr(i,j)+vr(i,j+1) ) * pi * dz * (r(j)+r(j+1))/2;
            end 
        end
    end
    
    % calculate Nnew -------------------------------------------------------------------
    
    % j=1
    for i = 2:size(N,1)
        Nnew(i,1) = N(i,1) - dt/V(i,1) * ( Fz(i,1) - Fz(i-1,1) + Fr(i,1) );
    end
    
    % j~=1  
    for i = 2:size(N,1)
        for j = 2:size(N,2)
            Nnew(i,j) =  N(i,j) - dt/V(i,j) * ( Fz(i,j) - Fz(i-1,j) + Fr(i,j) - Fr(i,j-1) );
        end
    end
    
    

    
    
    N = Nnew;  % update N
    t = t+dt;  % update t
    
    figure(2)
    imagesc(r,z,N)
    colorbar
    caxis([-8 8]*1e7)
    xlabel('r (m)')
    ylabel('z (m)')
    title(sprintf('charge density distribution at t = %.2f s',t));
    set(gca, 'Fontsize', 15)
    %pause
    
    F(n) = getframe(gcf) ;
    drawnow
    
    
    
end




% make movie ----------------------------------------------------------------

% create the video writer with 1 fps
writerObj = VideoWriter('myVideo.avi');
writerObj.FrameRate = 10;
  % set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
for i=1:length(F)
    % convert the image to a frame
    frame = F(i) ;    
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);






return

%% not used













Sz = zeros(1, size(ne,2)); 
Sr = zeros(1, size(ne,2)); 
for j = 2:size(ne,2)
    Sz(j) = pi*((r(j)+dr/2)^2-(r(j)-dr/2)^2);  % Sz2.5, Sz3.5...
    Sr(j) = 2*pi*(r(j)+dr/2)*dz;               % Sr2.5, Sr3.5...
end
Sz(1) = pi*(dr/2)^2;       % Sz1.5
Sr(1) = 2*pi*(dr/2)*dz;  % Sr1.5

Sz = repmat(Sz,size(ne,1),1);
Sr = repmat(Sr,size(ne,1),1);
V = Sz*dz;



% loop through time index
neNew = zeros(size(ne,1), size(ne,2));  % initialize it for every time step
for n = 1:ntimestep 

   
    % numerical solution


    for i = 1:size(ne,1)
       
        neNew(i,1) = ne(i,1) - dt/V(i,1)*( ne(i,1)*vr*Sr(i,1) );
        
%         uinp1 = u(i) - v*dt/dx * (u(i) - u(i-1));
%         uinp1Vec(i) = uinp1;

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
    xlabel('r')
    ylabel('z')
    title('electron density distribution')
    set(gca, 'Fontsize', 15)
    %pause(0.1)
    

%     
end





    
    
  

   


