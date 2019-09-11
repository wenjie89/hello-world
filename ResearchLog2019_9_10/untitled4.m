close all; clear; clc;
% one dimensional FCT algorithm with a Boris-Book limiter 


%% define parameters

% define spacial simulation domain -------------------------------------------------------
% dx = 0.05;

N = 50; % number of points
xmin  = 0;
xmax = 5;

dx = (xmax - xmin)/N;


x = xmin-dx : dx : xmax+dx; % introduce ghost nodes

% define parameters of time and velocity -------------------------------------------------

v = 1; % [m/s]
dt = 0.025; % [x]
t = 0; % current time [s]
T = 1; % maximum time [s]

% define the pulse type to be solved -----------------------------------------------------
u0 = rectPulse(x); % initial u0
figure(1)
plot(x,u0)
xlabel('x (m)')
ylabel('u(x,t)')
ylim([-0.5 1.5])
title('a rectangular pulse at t=0')
set(gca, 'Fontsize', 15)


%% the main algorithm

u = u0;

order = 6; % order of high order solution


ntimestep = floor(T/dt);


% loop through time index ----------------------------------------------------------------

for n = 1:ntimestep 
    
     % boundary conditions
    
 %   u(1) = u(3);
 %    u(2) = u(3);
 %   u(N+3) = u(N+1);
    
    
    delx = v*(n-1)*dt;  % spatial advancement corr. to time advancement
    
    % analytical solution
    ut  = rectPulse(x-delx); % transported pulse
    
    % numerical solution -----------------------------------------------------------------
    
    uinp1Vec = zeros(1,length(u)); % store the updated numerical solution 
    
   for i = 2:length(x)-1
        uinp1 = u(i) - v * dt/(2*dx) * (u(i+1)-u(i-1)) + 1/2 * (v * dt/dx)^2 * (u(i+1)-2*u(i)+u(i-1));  % the time advancement of u 
        uinp1Vec(i) = uinp1; % store it in a vector
    end
    
    % plot -------------------------------------------------------------------------------
    figure(2)
    plot(x,ut) % analytical solution
    hold on
    plot(x,u,'ro-')  % numerical solution
    xlabel('x (m)')
    ylabel('u(x,t)')
    ylim([-0.5 1.5])
    title(sprintf('FCT at t = %f',t))
    set(gca, 'Fontsize', 15)
  
    hold off
    pause(0.5)
   % saveas(gcf,sprintf('%d.png',n))

    % update solution --------------------------------------------------------------------
    u = uinp1Vec; % update the numerical solution in each loop 
    t = t+dt;     % update the current time 
    
end





    
    
  

   


