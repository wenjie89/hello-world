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
order = 6; % order of high order solution
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
    delx = v*(n-1)*dt;  % spatial advancement corr. to time advancement
    
    % analytical solution
    ut  = rectPulse(x-delx); % transported pulse
    
    % numerical solution
    uinp1Vec = zeros(1,length(x));
    
%     ================================================
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
            % second order flux
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
       % uitd   = u(i) - (FLiph - FLimh)/dx;                % eq. (6)       u_i^td
        uimotd = u(i-1) - (FLimh - FLimoh)/dx;             %               u_{i-1}^td
        uimttd = u(i-2) - (FLimoh - FLimth)/dx;            % 
        uipotd =   u(i+1) - (FLipoh - FLiph)/dx;           %               u_{i+1}^td
        uipttd =   u(i+2) - (FLipth - FLipoh)/dx;          %               u_{i+2}^td

       % step 5 -- Boris-Book limiter (might be wrong)
        if Aiph > 0 || Aiph == 0 
            Siph = 1;
        else 
            Siph = -1;
        end
        
        if Aimh > 0 || Aimh == 0 
            Simh = 1;
        else 
            Simh = -1;
        end
        
%         MatCPAiphC = [abs(Aiph), Siph*(uipttd-uipotd)*dx, Siph*(uitd-uimotd)*dx];
%         MatCPAimhC = [abs(Aimh), Simh*(uipotd-uitd)*dx, Simh*(uimotd-uimttd)*dx];

%         AiphC = Siph * max(0, min(MatCPAiphC) ); % eq. (7) (17)
%         AimhC = Simh * max(0, min(MatCPAimhC) );

       % step 6 -- limited antidiffusive fluxes
        uitd   = u(i) - ( 0.5*( v*(u(i+1)+u(i)) - abs(v)*(u(i+1)-u(i)) ) * dt -...
                          0.5*( v*(u(i)+u(i-1)) - abs(v)*(u(i)-u(i-1)) ) * dt )/dx;    
        %uinp1 = uitd-(AiphC-AimhC)/dx;             % the time advancement of u 
        uinp1 = uitd;   
        uinp1Vec(i) = uinp1;                       % store it in a vector
    
    end
%     ================================================


%     for i = 2:length(x)-1
%         uinp1 = u(i) - v * dt/(2*dx) * (u(i+1)-u(i-1)) + 1/2 * (v * dt/dx)^2 * (u(i+1)-2*u(i)+u(i-1));  % the time advancement of u 
%         uinp1Vec(i) = uinp1; % store it in a vector
%     end
%     
    figure(2)
    plot(x,ut) % analytical solution
    hold on
    plot(x,u,'ro-')  % numerical solution
    xlabel('x (m)')
    ylabel('u(x,t)')
    ylim([-0.5 1.5])
    title(sprintf('a rectangular pulse at t=%f',t))
    set(gca, 'Fontsize', 15)
    pause
    hold off
    %saveas(gcf,sprintf('LWcom at t = %.2f.png',t))

    
    u = uinp1Vec; % update the numerical solution in each loop 
    t = t+dt;     % update the current time index
    
end





    
    
  

   


