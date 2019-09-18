close all; clear; clc

choice = input('Enter a number: ');

% enter a number to choose the numerical scheme in the simulation

% 1. FOU
% 2. Donor cell
% 3. Lax-Fridrichs
% 4. Lax-Wendroff


%% define parameters

Nx = 50;  % number of grid points in x direction
Nt = 50;  % number of time step
CFLratio = 0.5;  % it cannot be greater than 1

v = 1;
xmin = 0;
xmax = 10;
dx = (xmax - xmin) / (Nx-1);  % cell size in x direction

dt = dx/v * CFLratio;
tmin = 0;
tmax = tmin + (Nt-1) * dt;
t = tmin:dt:tmax;


%% discretization of computational domain

x = xmin : dx : xmax;


u0 = rectPulse(x);  % initial pulst at t = 0;
figure(1)
plot(x,u0)

u = rectPulse(x);



%% transport the pulse

unp1 = u0; % wrong
unp1 = zeros(1,length(u)); % initialize u value at the next time index
totmass = zeros(1,Nt);  % initialize totmass function 

for n = 1:Nt  % loop through time
    
    % numerical algorithm
    for i = 2:length(u)  
        
        % calculate flux at interface   
        

    switch choice
        case 1
            % FOU-------------------------------------------
            Fiph = v*u(i);
            Fimh = v*u(i-1);
        case 2
            % Donor cell------------------------------------
            Fiph = v*u(i);
            Fimh = v*u(i-1);
        case 3
            disp('positive one')
    end

        % general expression of the FVM scheme
        unp1(i) = u(i) + dt/dx * ( Fimh - Fiph );  % update every spatial point of u
        
    end
    
    % analytical solution
    u_ana = rectPulse(x - v*(n-1)*dt);  % start from t = 0
    
    figure(2)
    plot(x,u,'bo-')
    hold on 
    plot(x,u_ana,'r-')
    hold off
    xlabel('x (m)')
    ylabel('u(x,t)')
    set(gca, 'FontSize', 15)
    if choice == 1
        legend('FOU', 'analytic');
    elseif choice == 2
         legend('Donor cell', 'analytic');
    elseif choice == 3
         legend('Lax-Fridrichs', 'analytic');
    else
        legend('Lax-Wendroff', 'analytic');
    end
   
    %pause(0.5)
    
    unp1mid = (unp1(1:length(unp1)-1) + unp1(2:length(unp1)))/2;  % mid-point value of unp1
    totmass(n) = sum(dx*unp1mid);   % calculate the total mass for each t index
    
    
    % t = tmin + n * dt;  % track current time
    u = unp1; % update u 
    
end


figure(3)
clf
plot(t(1:end),totmass(1:end),'x')
ylim([0 1.2*max(totmass)])
xlabel('t (s)')
ylabel('total mass')
set(gca,'Fontsize', 15)



















