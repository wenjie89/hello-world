close all; clear; clc;

R =[0.3e-3 0.5e-3];
E0 = 4e5; 
eps0 = 8.85e-12;
Qdr = (0:100:500)*1e-12;
%z = linspace(1e-9,1e-1,100);
z = 10.^(-7:0.1:-1);

Ecell = cell(2);
for j = 1:length(R)
    E = zeros(length(Qdr),length(z));
    for i = 1:length(Qdr)
        E(i,:) = E0*(1+2*R(j)^3./(R(j)+z).^3) + Qdr(i)./(4*pi*eps0*(R(j)+z).^2);
    end
    Ecell{j} = E;
end

E1 = Ecell{1};
E2 = Ecell{2};

%%

figure(1)  % electric field as a function of z for R = 0.25 mm
clf
semilogx(z*1e3,E1(1,:)*1e-5,'k')
hold on 
semilogx(z*1e3,E1(2,:)*1e-5,'k')
semilogx(z*1e3,E1(3,:)*1e-5,'k')
semilogx(z*1e3,E1(4,:)*1e-5,'k')
semilogx(z*1e3,E1(5,:)*1e-5,'k')
semilogx(z*1e3,E1(6,:)*1e-5,'k')

%legend('Q_{dr} = 0','Q_{dr} = 100 pC', 'Q_{dr} = 200 pC', 'Q_{dr} = 300 pC',...
%    'Q_{dr} = 400 pC', 'Q_{dr} = 500 pC')
xlabel('z (mm)')
ylabel('Electric field (kV/cm)')
set(gca,'FontSize',15)

figure(2)  % electric field as a function of z for R = 0.75 mm
clf
semilogx(z*1e3,E2(1,:)*1e-5,'k')
hold on 
semilogx(z*1e3,E2(2,:)*1e-5,'k')
semilogx(z*1e3,E2(3,:)*1e-5,'k')
semilogx(z*1e3,E2(4,:)*1e-5,'k')
semilogx(z*1e3,E2(5,:)*1e-5,'k')
semilogx(z*1e3,E2(6,:)*1e-5,'k')

%legend('Q_{dr} = 0','Q_{dr} = 100 pC', 'Q_{dr} = 200 pC', 'Q_{dr} = 300 pC',...
%    'Q_{dr} = 400 pC', 'Q_{dr} = 500 pC')
xlabel('z (mm)')
ylabel('Electric field (kV/cm)')
set(gca,'FontSize',15)

%%
alpha_air1_E1 = zeros(size(E1)); % effective ionization coef. for R = 0.25mm
alpha_morrow_E1 = zeros(size(E1));
alpha_air1_E2 = zeros(size(E2)); % effective ionization coef. for R = 0.75mm
alpha_morrow_E2 = zeros(size(E2));

for i = 1:size(E1,1)
    nui_air1=air1(E1(i,:),0,10);   % ionization freq
    nua_air1=air1(E1(i,:),0,2);    % attachment freq
    mue_air1 = air1(E1(i,:),0,11); % electron mobility
    ve_air1 = -mue_air1.*E1(i,:);  % electron drift velocity
    alpha_i_air1 = nui_air1./abs(ve_air1); % ionization coefficient
    alpha_a_air1 = nua_air1./abs(ve_air1); % attachment coefficient
    alpha_air1 = alpha_i_air1 - alpha_a_air1;   % effective ionization coef.
    alpha_air1_E1(i,:) = alpha_air1;

    nui_morrow=morrowair(E1(i,:),0,1);    % ionization freq
    nua_morrow=morrowair(E1(i,:),0,2);    % attachment freq
    mue_morrow = morrowair(E1(i,:),0,4);  % electron mobility
    ve_morrow = -mue_morrow.*E1(i,:);     % electron drift velocity
    alpha_i_morrow = nui_morrow./abs(ve_morrow);  % ionization coefficient
    alpha_a_morrow = nua_morrow./abs(ve_morrow);  % attachment coefficient
    alpha_morrow = alpha_i_morrow - alpha_a_morrow;  % effective ionization coef.
    alpha_morrow_E1(i,:) = alpha_morrow; 
end

figure(3)
clf
for i = 1:size(E1,1)
semilogx(z*1e3,alpha_air1_E1(i,:),'k')
hold on
end
xlabel('z (mm)')
ylabel(' \alpha_{ion} - \alpha_{att} (m^{-1})')
%legend('Q_{dr} = 0','Q_{dr} = 100 pC', 'Q_{dr} = 200 pC', 'Q_{dr} = 300 pC',...
%    'Q_{dr} = 400 pC', 'Q_{dr} = 500 pC')
set(gca,'FontSize',15)

figure(4)
clf
for i = 1:size(E1,1)
semilogx(z*1e3,alpha_morrow_E1(i,:),'k')
hold on
end
xlabel('z (mm)')
ylabel(' \alpha_{ion} - \alpha_{att} (m^{-1})')
%legend('Q_{dr} = 0','Q_{dr} = 100 pC', 'Q_{dr} = 200 pC', 'Q_{dr} = 300 pC',...
%    'Q_{dr} = 400 pC', 'Q_{dr} = 500 pC')
set(gca,'FontSize',15)


for i = 1:size(E2,1)
    nui_air1=air1(E2(i,:),0,10);   % ionization freq
    nua_air1=air1(E2(i,:),0,2);    % attachment freq
    mue_air1 = air1(E2(i,:),0,11); % electron mobility
    ve_air1 = -mue_air1.*E2(i,:);  % electron drift velocity
    alpha_i_air1 = nui_air1./abs(ve_air1); % ionization coefficient
    alpha_a_air1 = nua_air1./abs(ve_air1); % attachment coefficient
    alpha_air1 = alpha_i_air1 - alpha_a_air1;   % effective ionization coef.
    alpha_air1_E2(i,:) = alpha_air1;

    nui_morrow=morrowair(E2(i,:),0,1);    % ionization freq
    nua_morrow=morrowair(E2(i,:),0,2);    % attachment freq
    mue_morrow = morrowair(E2(i,:),0,4);  % electron mobility
    ve_morrow = -mue_morrow.*E2(i,:);     % electron drift velocity
    alpha_i_morrow = nui_morrow./abs(ve_morrow);  % ionization coefficient
    alpha_a_morrow = nua_morrow./abs(ve_morrow);  % attachment coefficient
    alpha_morrow = alpha_i_morrow - alpha_a_morrow;  % effective ionization coef.
    alpha_morrow_E2(i,:) = alpha_morrow; 
end

figure(5)
clf
for i = 1:size(E2,1)
semilogx(z*1e3,alpha_air1_E2(i,:),'k')
hold on
end
xlabel('z (mm)')
ylabel(' \alpha_{ion} - \alpha_{att} (m^{-1})')
%legend('Q_{dr} = 0','Q_{dr} = 100 pC', 'Q_{dr} = 200 pC', 'Q_{dr} = 300 pC',...
%    'Q_{dr} = 400 pC', 'Q_{dr} = 500 pC')
ylim([0 3.1e5])
xlim([1e-4 1e2])
set(gca,'FontSize',15)

figure(6)
clf
for i = 1:size(E2,1)
semilogx(z*1e3,alpha_morrow_E2(i,:),'k')
hold on
end
xlabel('z (mm)')
ylabel(' \alpha_{ion} - \alpha_{att} (m^{-1})')
%legend('Q_{dr} = 0','Q_{dr} = 100 pC', 'Q_{dr} = 200 pC', 'Q_{dr} = 300 pC',...
%    'Q_{dr} = 400 pC', 'Q_{dr} = 500 pC')
set(gca,'FontSize',15)

return
%% not used

figure(7)
clf
plot(z,alpha_air1_E1(1,:))
hold on
plot(z,alpha_air1_E1(2,:))
plot(z,alpha_air1_E1(3,:))
plot(z,alpha_air1_E1(4,:))
plot(z,alpha_air1_E1(5,:))
plot(z,alpha_air1_E1(6,:))
xlim([0 1e-3])

alpha_air1_E1(alpha_air1_E1<0)=0; % remove the ones with negative alpha

sum(alpha_air1_E1(2,:)*(z(2)-z(1)))








figure(8)
clf
plot(z,alpha_morrow_E1(1,:))
hold on
plot(z,alpha_morrow_E1(2,:))
plot(z,alpha_morrow_E1(3,:))
plot(z,alpha_morrow_E1(4,:))
plot(z,alpha_morrow_E1(5,:))
plot(z,alpha_morrow_E1(6,:))
xlim([0 1e-3])

figure(9)
clf
plot(z,alpha_air1_E2(1,:))
hold on
plot(z,alpha_air1_E2(2,:))
plot(z,alpha_air1_E2(3,:))
plot(z,alpha_air1_E2(4,:))
plot(z,alpha_air1_E2(5,:))
plot(z,alpha_air1_E2(6,:))
xlim([0 1e-3])


figure(10)
clf
plot(z,alpha_morrow_E2(1,:))
hold on
plot(z,alpha_morrow_E2(2,:))
plot(z,alpha_morrow_E2(3,:))
plot(z,alpha_morrow_E2(4,:))
plot(z,alpha_morrow_E2(5,:))
plot(z,alpha_morrow_E2(6,:))
xlim([0 1e-3])













