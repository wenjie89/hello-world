close all; clear; clc;
% one dimensional FCT algorithm with a Boris-Book limiter 
% (1) Initialized uitd etc. to multidimension.
% (2) Donor cell test is successful.
% (3) Conservative test is successful.

order = input('Enter the order (2, 4 or 6) of the high order algorithm: '); % order of high order solution

case_N = input('Choose the low order scheme (1 for eq.(14), 2 for donor cell): ');

pulse_N = input('Choose the pulse type (1 for rectangular pulse, 2 for Gaussian pulse): ');

test_N = input('Only test the low order shceme (1 for the whole FCT, 2 only for the low order scheme, 3 only for the high order scheme) :');

% max_min_N = input('Choose the max and min value for the Zalesak limiter (1 for a safe choice, 2 for a better choice):');



%% define parameters

% define spacial simulation domain -------------------------------------------------------

N = 100; % number of points
xmin  = 0;
xmax = 3;

dx = (xmax - xmin)/N;


%x = xmin-dx : dx : xmax+dx; % introduce ghost nodes
x = xmin : dx : xmax;

% define parameters of time and velocity -------------------------------------------------

CFL = 0.13;
dt = 0.015; % [x]
v = CFL*dx/dt; % [m/s]
t = 0; % current time [s]
T = 1.8-dt; % maximum time [s]
tVec = t:dt:T;


% define the pulse type to be solved -----------------------------------------------------
switch pulse_N
    case 1  % rectangular pulse
        u0 = rectPulse(x); % initial u0
    case 2
        u0 = exp(-200*(x-0.25).^2); % Gaussian pulse 
end
        
    


figure(1)
plot(x,u0)
xlabel('x (m)')
ylabel('u(x,t)')
xlim([0 3])
% ylim([-0.5 1.5])
title('a rectangular pulse at t=0')
set(gca, 'Fontsize', 15)


%% the main algorithm

u = u0;


ntimestep = floor(T/dt);


% initialization
uitd   = zeros(1,length(u));
uimotd = zeros(1,length(u));
uimttd = zeros(1,length(u));
uipotd = zeros(1,length(u));
uipttd = zeros(1,length(u));
uinp1  = zeros(1,length(u));
uiHtd = zeros(1,length(u));  % store high order solution
totmass = zeros(1,ntimestep);



% loop through time index ----------------------------------------------------------------

for n = 1:ntimestep 
    
     % boundary conditions
    
 %   u(1) = u(3);
 %    u(2) = u(3);
 %   u(N+3) = u(N+1);
    
    
    delx = v*n*dt;  % spatial advancement corr. to time advancement
    
    % analytical solution
    switch pulse_N
        case 1  % rectangular pulse 
            ut  = rectPulse(x-delx); % transported pulse
        case 2  % Gaussian pulse 
            ut = exp(-200*(x-delx-0.25).^2);
    end
    
    
    % numerical solution -----------------------------------------------------------------
    
    % uinp1Vec = zeros(1,length(u)); % store the updated numerical solution 
    
    % for i = 4:length(u)-3
    for i = 3:length(u)-2    
        
        % main algorithm of FCT ----------------------------------------------------------
        
        
        % step 1 -- low order scheme

        switch case_N
            case 1    % eq. (14)
                FLiph =  0.5*( v*(u(i+1)+u(i)) - abs(v)*(u(i+1)-u(i)) ) * dt;   
                FLimh =  0.5*( v*(u(i)+u(i-1)) - abs(v)*(u(i)-u(i-1)) ) * dt; 
                FLimoh = 0.5*( v*(u(i-1)+u(i-2)) - abs(v)*(u(i-1)-u(i-2)) ) * dt;
               % FLimth = 0.5*( v*(u(i-2)+u(i-3)) - abs(v)*(u(i-2)-u(i-3)) ) * dt;
                FLipoh = 0.5*( v*(u(i+2)+u(i+1)) - abs(v)*(u(i+2)-u(i+1)) ) * dt; 
               % FLipth = 0.5*( v*(u(i+3)+u(i+2)) - abs(v)*(u(i+3)-u(i+2)) ) * dt; 
            case 2   % donor cell
                FLiph =  v*u(i) * dt;
                FLimh =  v*u(i-1) * dt;
                FLimoh = v*u(i-2) * dt;
               % FLimth = v*u(i-3) * dt;
                FLipoh = v*u(i+1) * dt;
               % FLipth = v*u(i+2) * dt;             
        end

        % step 2 -- high order scheme
        if order == 2
            %  second order flux
            FHiph =   0.5*v*(u(i+1)+u(i)) * dt;  
            FHimh =   0.5*v*(u(i)+u(i-1)) * dt;
            FHipoh =  0.5*v*(u(i+2)+u(i+1)) * dt;
            FHimoh =  0.5*v*(u(i-1)+u(i-2)) * dt;
        elseif order == 4
             %  forth order flux eq. (15)
            FHiph =  ( 7/12*v*(u(i+1)+u(i)) - 1/12*v*(u(i+2)+u(i-1)) ) * dt;  
            FHimh =  ( 7/12*v*(u(i)+u(i-1)) - 1/12*v*(u(i+1)+u(i-2)) ) * dt;
            FHipoh = ( 7/12*v*(u(i+2)+u(i+1)) - 1/12*v*(u(i+3)+u(i)) ) * dt;
            FHimoh = ( 7/12*v*(u(i-1)+u(i-2)) - 1/12*v*(u(i)+u(i-3)) ) * dt;
        elseif order == 6
             %  sixth order flux
            FHiph  = ( 37/60*v*(u(i+1)+u(i)) - 2/15*v*(u(i+2)+u(i-1)) +1/60*v*(u(i+3)+u(i-2)) ) * dt; 
            FHimh  = ( 37/60*v*(u(i)+u(i-1)) - 2/15*v*(u(i+1)+u(i-2)) +1/60*v*(u(i+2)+u(i-3)) ) * dt;  
            FHipoh = ( 37/60*v*(u(i+2)+u(i+1)) - 2/15*v*(u(i+3)+u(i)) +1/60*v*(u(i+4)+u(i-1)) ) * dt; 
            FHimoh = ( 37/60*v*(u(i-1)+u(i-2)) - 2/15*v*(u(i)+u(i-3)) +1/60*v*(u(i+1)+u(i-4)) ) * dt;  
        else
            error('The order is wrong!')
        end

        % step 3 -- antidiffusion
        Aiph   = FHiph - FLiph;  % eq. (5)
        Aimh   = FHimh - FLimh;
        Aipoh  = FHipoh - FLipoh;
        Aimoh  = FHimoh - FLimoh;

        % step 4 -- compute updated low order solution
        uitd(i)   = u(i) - (FLiph - FLimh)/dx;                % eq. (6)       u_i^td
%         uimotd(i) = u(i-1) - (FLimh - FLimoh)/dx;             %               u_{i-1}^td
%         uimttd(i) = u(i-2) - (FLimoh - FLimth)/dx;            % 
%         uipotd(i) =   u(i+1) - (FLipoh - FLiph)/dx;           %               u_{i+1}^td
%         uipttd(i) =   u(i+2) - (FLipth - FLipoh)/dx;          %               u_{i+2}^td
        uiHtd(i) = u(i) - (FHiph - FHimh)/dx; 

        % step 5 -- Boris-Book limiter (might be wrong)
%         if Aiph > 0 || Aiph == 0 
%             Siph = 1;
%         else 
%             Siph = -1;
%         end
%         
%         if Aimh > 0 || Aimh == 0 
%             Simh = 1;
%         else 
%             Simh = -1;
%         end
%         
%         MatCPAiphC = [abs(Aiph), Siph*(uipttd(i)-uipotd(i))*dx, Siph*(uitd(i)-uimotd(i))*dx];
%         MatCPAimhC = [abs(Aimh), Simh*(uipotd(i)-uitd(i))*dx, Simh*(uimotd(i)-uimttd(i))*dx];
% 
%         AiphC = Siph * max(0, min(MatCPAiphC) ); % eq. (7) (17)
%         AimhC = Simh * max(0, min(MatCPAimhC) );
        
        % Zalesak limiter--------------------------------------------------
        
      
        % (1) select uimax and uimin
        
        % choice one
%         uimax =  max([uimotd(i), uitd(i), uipotd(i)]);
%         uipomax = max([uimotd(i+1), uitd(i+1), uipotd(i+1)]);
%         uimin = min([uimotd(i), uitd(i), uipotd(i)]);
%         uipomin = min([uimotd(i+1), uitd(i+1), uipotd(i+1)]);
        
        % a better choice
        uia =   max([u(i), uitd(i)]);
        uimoa = max([u(i-1), uitd(i-1)]);
        uipoa = max([u(i+1), uitd(i+1)]);
        uipta = max([u(i+2), uitd(i+2)]);
        
        uimax = max([uimoa,uia,uipoa]);
        uipomax = max([uia,uipoa,uipta]);
        
        uib =   min([u(i), uitd(i)]);
        uimob = min([u(i-1), uitd(i-1)]);
        uipob = min([u(i+1), uitd(i+1)]);
        
        
        uimin = min([uimoa,uia,uipoa]);
        uipomin = min([uia,uipoa,uipta]);
        
        
        
        % (2) build P, Q, and R
        P_plus(i) = max(0,Aimh) - min(0,Aiph);
        Q_plus(i) = (uimax - uitd(i))*dx;
        if P_plus(i) > 0
            R_plus(i) = min(1,Q_plus(i)/P_plus(i));
        elseif P_plus(i) == 0
                 R_plus(i) = 0;
        else
            error('Value of P_plus(i) is wrong!')
        end
        
        
        P_plus(i+1) = max(0,Aiph) - min(0,Aipoh);
        Q_plus(i+1) = (uipomax - uitd(i+1))*dx;
        if P_plus(i+1) > 0
            R_plus(i+1) = min(1,Q_plus(i+1)/P_plus(i+1));
        elseif P_plus(i+1) == 0
                 R_plus(i+1) = 0;
        else
            error('Value of P_plus(i+1) is wrong!')
        end
        
        
        P_minus(i) = max(0,Aiph) - min(0,Aimh);
        Q_minus(i) = (uitd(i) - uimin)*dx;
        if P_minus(i) > 0
            R_minus(i) = min(1,Q_minus(i)/P_minus(i));
        elseif P_minus(i) == 0
            R_minus(i) = 0;
        else 
            error('Value of P_minus(i) is wrong!')
        end
        
        
        P_minus(i+1) = max(0,Aipoh) - min(0,Aimoh);
        Q_minus(i+1) = (uitd(i+1) - uipomin)*dx;
        if P_minus(i+1) > 0
            R_minus(i+1) = min(1,Q_minus(i+1)/P_minus(i+1));
        elseif P_minus(i+1) == 0
            R_minus(i+1) = 0;
        else 
            error('Value of P_minus(i+1) is wrong!')
        end
        
        
        
        % (3) compute the coefficient of the limiter 
        if Aiph > 0 || abs(Aiph/FLiph) < 1e-10
            Ciph = min(R_plus(i+1),R_minus(i));
        else
            Ciph = min(R_plus(i),R_minus(i+1));
        end
        
        if Aimh > 0 || abs(Aimh/FLimh) < 1e-10
            Cimh = min(R_plus(i),R_minus(i-1));
        else
            Cimh = min(R_plus(i-1),R_minus(i));
        end
              
        
        if Ciph > 1+1e-10 || Ciph < -1e-10 ||  Cimh > 1+1e-10 || Cimh < -1e-10  % check right / wrong
            error('Wrong value for Ciph!')
        end
        
        % (4) compute AiphC and AimhC
        
        AiphC = Ciph*Aiph;
        AimhC = Cimh*Aimh;
        

        % step 6 -- limited antidiffusive fluxes
        
        switch test_N 
            case 1   % the whole FCT scheme
                uinp1(i) = uitd(i)-(AiphC-AimhC)/dx;   % the time advancement of u 
            case 2   % only the low order scheme
                uinp1(i) = uitd(i);   
            case 3   % only the high order scheme
                uinp1(i) = uiHtd(i); 
        end      
        
        % uinp1Vec(i) = uinp1;                           % store it in a vector
    
    end
    
    % plot -------------------------------------------------------------------------------
    figure(2)
    plot(x,ut) % analytical solution
    hold on
    plot(x(1:end),u(1:end),'ro-')  % numerical solution
    xlabel('x (m)')
    ylabel('u(x,t)')
    ylim([-0.5 1.5])
    title(sprintf('FCT at t = %f',t))
    set(gca, 'Fontsize', 15)
    hold off
    pause(0.1)
   % saveas(gcf,sprintf('%d.png',n))
   
    % conservative check
    uinp1mid = (uinp1(1:length(uinp1)-1) + uinp1(2:length(uinp1)))/2;  % mid-point value of unp1
    totmass(n) = sum(dx*uinp1mid);   % calculate the total mass for each t index

    % update solution --------------------------------------------------------------------
    u = uinp1;    % update the numerical solution in each loop 
    t = t+dt;     % update the current time 
    
end

figure(3)
clf
plot(tVec(1:end-1),totmass,'x')
ylim([0 1.2*max(totmass)])
xlabel('t (s)')
ylabel('total mass')
set(gca,'Fontsize', 15)




    
    
  

   


