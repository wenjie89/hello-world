close all; clear; clc;
% one dimensional FCT algorithm with a Boris-Book limiter 


%% define parameters

% define spacial simulation domain -------------------------------------------------------
% dx = 0.05;

N = 100; % number of points
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
    u(N+3) = u(N+1);
    
    
    delx = v*n*dt;  % spatial advancement corr. to time advancement
    
    % analytical solution
    ut  = rectPulse(x-delx); % transported pulse
    
    % numerical solution -----------------------------------------------------------------
    
    uinp1Vec = zeros(1,length(u)); % store the updated numerical solution 
    
    for i = 4:length(u)-3
        
        % main algorithm of FCT ----------------------------------------------------------
        
        % step 1 -- low order scheme
        FLiph =  0.5*( v*(u(i+1)+u(i)) - abs(v)*(u(i+1)-u(i)) ) * dt;    % eq. (14)
        FLimh =  0.5*( v*(u(i)+u(i-1)) - abs(v)*(u(i)-u(i-1)) ) * dt; 
        FLimoh = 0.5*( v*(u(i-1)+u(i-2)) - abs(v)*(u(i-1)-u(i-2)) ) * dt;
        FLimth = 0.5*( v*(u(i-2)+u(i-3)) - abs(v)*(u(i-2)-u(i-3)) ) * dt;
        FLipoh = 0.5*( v*(u(i+2)+u(i+1)) - abs(v)*(u(i+2)-u(i+1)) ) * dt; 
        FLipth = 0.5*( v*(u(i+3)+u(i+2)) - abs(v)*(u(i+3)-u(i+2)) ) * dt; 

        % step 2 -- high order scheme
        if order == 2
            %  second order flux
            FHiph =  0.5*v*(u(i+1)+u(i)) * dt;  
            FHimh =  0.5*v*(u(i)+u(i-1)) * dt;
        elseif order == 4
             %  forth order flux eq. (15)
            FHiph = ( 7/12*v*(u(i+1)+u(i)) - 1/12*v*(u(i+2)+u(i-1)) ) * dt;  
            FHimh = ( 7/12*v*(u(i)+u(i-1)) - 1/12*v*(u(i+1)+u(i-2)) ) * dt;
        elseif order == 6
             %  sixth order flux
            FHiph = ( 37/60*v*(u(i+1)+u(i)) - 2/15*v*(u(i+2)+u(i-1)) +1/60*v*(u(i+3)+u(i-2)) ) * dt; 
            FHimh = ( 37/60*v*(u(i)+u(i-1)) - 2/15*v*(u(i+1)+u(i-2)) +1/60*v*(u(i+2)+u(i-3)) ) * dt;  
        else
            error('The order is wrong!')
        end

        % step 3 -- antidiffusion
        Aiph   = FHiph - FLiph;  % eq. (5)
        Aimh   = FHimh - FLimh;

        % step 4 -- compute updated low order solution
        uitd   = u(i) - (FLiph - FLimh)/dx;                % eq. (6)       u_i^td
        uimotd = u(i-1) - (FLimh - FLimoh)/dx;             %               u_{i-1}^td
        uimttd = u(i-2) - (FLimoh - FLimth)/dx;            % 
        uipotd =   u(i+1) - (FLipoh - FLiph)/dx;           %               u_{i+1}^td
        uipttd =   u(i+2) - (FLipth - FLipoh)/dx;          %               u_{i+2}^td

        % step 5 -- Boris-Book limiter 
        if Aiph > 0 || Aiph == 0 
            S = 1;
        else 
            S = -1;
        end
        
        MatCPAiphC = [abs(Aiph), S*(uipttd-uipotd)*dx, S*(uitd-uimotd)*dx];
        MatCPAimhC = [abs(Aimh), S*(uipotd-uitd)*dx, S*(uimotd-uimttd)*dx];

        AiphC = S * max(0, min(MatCPAiphC) ); % eq. (7) (17)
        AimhC = S * max(0, min(MatCPAimhC) );

        % step 6 -- limited antidiffusive fluxes
        uinp1 = uitd-(AiphC-AimhC)/dx;             % the time advancement of u 
        uinp1Vec(i) = uinp1;                       % store it in a vector
    
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
    pause
    hold off
   % saveas(gcf,sprintf('%d.png',n))

    % update solution --------------------------------------------------------------------
    u = uinp1Vec; % update the numerical solution in each loop 
    t = t+dt;     % update the current time 
    
end





    
    
  

   


