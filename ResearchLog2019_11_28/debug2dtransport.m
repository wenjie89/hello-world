%% 2d donor cell algorighm for the transport equation

close all; clear; clc;
delete *.avi

% define parameters--------------------------------------------------------

eps0 = 8.85e-12;    % vacuum permittivity
rou0 = 20e-12;    % charge density in the sphere [C/m^3]  
q = 1.602e-19;      % elementary charge [C] 

Lz = 9e-2;          % Lx = 3cm
Lr = 3e-2;          % Ly = 3cm

z0 = Lz/2;          % z coordinate of the center of the source sphere
r0 = Lr/6;          % radius of the source sphere

nr = 201;           % number of points 

w = 1.9;            % relaxation factor for SOR
epsilon = 1e-8;     % accuracy factor for SOR

h = 0;              % altitude

% E_add = 2e-7;       % environment electric field [V/m]

E_add = 0;



% discretization domain----------------------------------------------------

dr = Lr/(nr-1);
dz = dr;

% z = 0:dz:Lz;
% r = 0:dr:Lr;

z = dz/2:dz:Lz; % z1, z2, z3...
r = 0:dr:Lr;    % r1, r2, r3...


[R, Z] = meshgrid(r,z);




%% numerical solution


% define the source sphere-------------------------------------------------
S = rou0 * ones(size(Z,1),size(Z,2));  %  [C/m^3]

sourceRegion = (Z-z0).^2 + R.^2 - r0^2; 

[idxX, idxY] = find(sourceRegion > 0);

for i = 1:length(idxX)
    S(idxX(i),idxY(i)) = 0;
end

N = S/q;   % [m^-3] density distribution of positive charges


% volume elements for the whole computational domain ----------------------
V = zeros(size(Z,1),size(Z,2));  
V(:,1) = pi*(dr/2)^2*dz;
for j=2:size(V,2)
    V(:,j) = pi * ( (r(j)+dr/2 )^2 - ( r(j-1)+dr/2 )^2 )  * dz;
end


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

% define parameters -------------------------------------------------------

% CFLz = 0.35;
% CFLz = 0.35;
% CFLr = 0.16;
% CFLr = 0;



% vz = CFLz*dz/dt;  % [m/s]
% vr = CFLr*dr/dt; 

dt = 0.02; % [s]
t = 0; % current time [s]
T = 20; % maximum time [s]



% transport the pulse

ntimestep = floor(T/dt);

% check conservation
TotCharge = zeros(1,ntimestep+1);
TotCharge(1) = sum(sum(N.*V));



% z = dz/2:dz:Lz; % z1, z2, z3...
% r = 0:dr:Lr;    % r1, r2, r3...


% % volume elements of the whole computational domain -----------------------
% 
% V  = zeros(size(N,1),size(N,2));                        % initialize V_ij
% V(:,1) = pi * ((r(1)+r(2))/2)^2 * dz;                   % V_i1
% for j = 2:size(V,2)
%     V(:,j) = pi * ((r(j)+dr/2)^2-(r(j)-dr/2)^2) * dz;   % V_ij
% end

% figure
% plot(V(1,:))




% % assume constant velocity in both directions
% 
% vr = CFLr*dr/dt * ones(size(N,1),size(N,2)); 
% vz = CFLz*dz/dt * ones(size(N,1),size(N,2)); 


Fz = zeros(size(N,1),size(N,2));
Fr = zeros(size(N,1),size(N,2));


% Nnew = zeros(size(N,1),size(N,2));  % initialize Nnew matrix
Nnew = N;  % initialize Nnew matrix
Snew = S;  % initialize Snew matrix

timeVec = linspace(0,T,ntimestep+1);
positiveCheck = zeros(1,length(timeVec));
itVec = zeros(1,length(timeVec)-1);
afVec = zeros(1,length(timeVec)-1);
%potCell = cell(1,length(timeVec)-1);  % potential cell
vrCell = cell(1,length(timeVec)-1);
vzCell = cell(1,length(timeVec)-1);





%      figure
%      quiver(R(1:10:end,1:10:end),Z(1:10:end,1:10:end),vr(1:10:end,1:10:end),vz(1:10:end,1:10:end))
    
    
%     figure(2)
%     imagesc(r,z,v)
%     colorbar
%     
%     
%     figure(3)
%     imagesc(r,z,Er)
%     colorbar
%     
%     
%     figure(4)
%     imagesc(r,z,vr)
%     colorbar
    
    

 

tic
for n = 1:ntimestep
    
    % calculate velocity --------------------------------------------------
    
%     [v, it] = CalcPotential(r, z, Lr, Lz, S, V, w, epsilon);
%     itVec(n) = it;
%     [Er, Ez] = gradient(v); 
%     Er = -Er/dr;
%     Ez = -Ez/dz;
    
    
    [v, it, af] = CalcPotential(r, z, Lr, Lz, S, V, w, epsilon);
    %potCell{n} = v;
    itVec(n) = it;
    afVec(n) = af;
    [Er, Ez] = gradient(v); 
    Er = -Er/dr;
    Ez = -Ez/dz+E_add;
    

    vr = morrowair(Er,h,4) .* Er;     % drift velocity of positive charges
    vz = morrowair(Ez,h,4) .* Ez;
    
    vrCell{n} = vr;
    vzCell{n} = vz;
    
    
    fprintf('Velocity calc. at time %f is done.\n', t);
    fprintf('CFLr = %.3f.\n', max(max(vr))*dt/dr);
    fprintf('CFLz = %.3f.\n', max(max(vz))*dt/dz);
    fprintf('maximal vr = %.3f.\n',max(max(vr)));
    fprintf('maximal vz = %.3f.\n',max(max(vz)));
    
    
    
    
    
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
    
    
    % update variables -----------------------------------------------------------------
    
    N = Nnew;    % update N
    S = Nnew*q;  % update S
    t = t+dt;    % update t  
    
    
    % positive check -------------------------------------------------------------------
    
    if  ~isempty(find(N<0, 1)) 
        positiveCheck(n+1)=1;
        warning('Negative density appears!')
    end
   

    % calculate total charge ------------------------------------------------------------
     
    TotCharge(n+1) = sum(sum(N.*V));
    
    % make plots ------------------------------------------------------------------------
    
    figure(2)
    imagesc(r,z,N)
    colorbar
    caxis([-8 8]*1e7)
    xlabel('r (m)')
    ylabel('z (m)')
    title(sprintf('charge density distribution at t = %.2f s',t));
    set(gca, 'Fontsize', 15)
    shg
    
%     figure(3)
%     quiver(R(1:50:end),Z(1:50:end),vr(1:50:end),vz(1:50:end))
%     xlabel('vr (m/s)')
%     ylabel('vz (m/s)')
%     title(sprintf('velocity n at t = %.2f s',t));
%     shg
    

    
    F(n) = getframe(gcf) ;
    drawnow
    
    
    
end

toc



% figure
% plot(timeVec,TotCharge)
% xlim([0 5])


if ~isempty(find(positiveCheck~=0, 1))
    warning('Positive check is fails.')
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








