close all; clear; clc;
% It works well although it may not be totally correct.

%% define a rectangular pulse at t = 0

dx = 0.05;
xmin  = 0;
xmax = 5;
x = xmin:dx:xmax;
u = rectPulse(x); % initial u
figure(1)
plot(x,u)
xlabel('x (m)')
ylabel('u(x,t)')
ylim([-0.5 1.5])
title('a rectangular pulse at t=0')
set(gca, 'Fontsize', 15)

%% the main algorithm

% define parameters --------------------------------------------------

v = 1; % [m/s]
dt = 0.025; % [x]
% dx = 1e-2; % [m]
t = 0; % current time [s]
T = 1; % maximum time [s]
N = length(x); % number of points


% transport the pulse

ntimestep = T/dt;

% loop through time index
for n = 1:ntimestep 
    delx = v*n*dt;  % spatial advancement corr. to time advancement
    
    % analytical solution
    ut  = rectPulse(x-delx); % transported pulse
    
    % numerical solution
    uinp1Vec = zeros(1,length(x));
    for i = 2:length(x)-1
        uinp1 = u(i) - v * dt/(2*dx) * (u(i+1)-u(i-1)) + 1/2 * (v * dt/dx)^2 * (u(i+1)-2*u(i)+u(i-1));  % the time advancement of u 
        uinp1Vec(i) = uinp1; % store it in a vector
    end
    
    figure(2)
    plot(x,ut) % analytical solution
    hold on
    plot(x,u,'ro-')  % numerical solution
    xlabel('x (m)')
    ylabel('u(x,t)')
    ylim([-0.5 1.5])
    title(sprintf('a rectangular pulse at t=%f',t))
    set(gca, 'Fontsize', 15)
  %  pause(0.5)
    hold off
    %saveas(gcf,sprintf('LWcom at t = %.2f.png',t))

    
    u = uinp1Vec; % update the numerical solution in each loop 
    t = t+dt;     % update the current time index
    
end





    
    
  

   


