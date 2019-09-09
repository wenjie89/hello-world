close all; clear all; clc;

e=10.^(5:0.00001:10);


h = fzero( @(h) satm(h)-0.4, 8); % find the altitude corresponding to 0.4 atm
% h = 0;

% calculate nu first
nui_air1 = air1(e,h,10);
nui_morrow = morrowair(e,h,1);


nua_air1 = air1(e,h,2);
nua_morrow = morrowair(e,h,2);


nu_air1 = nui_air1 - nua_air1;
nu_morrow = nui_morrow - nua_morrow;

mue_air1 = air1(e,h,11);
mue_morrow = morrowair(e,h,4);

% % calculate the field at which nu_air1 = nu_morrow
% fzero(@(x) (air1(x,h,10)- air1(x,h,2)) - (morrowair(x,h,1)- morrowair(x,h,2)), 1e7 )



% calculate the drift velocity then
ve_air1 = -mue_air1.*e;
ve_morrow = -mue_morrow.*e;


% obtain alpha from nu and ve
alpha_i_air1 = nui_air1./abs(ve_air1);
alpha_a_air1 = nua_air1./abs(ve_air1);

alpha_i_morrow = nui_morrow./abs(ve_morrow);
alpha_a_morrow = nua_morrow./abs(ve_morrow);

alpha_air1 = alpha_i_air1- alpha_a_air1;
alpha_morrow = alpha_i_morrow - alpha_a_morrow;

E_eq = fzero(@(x) (air1(x,h,10) - air1(x,h,2)) - (morrowair(x,h,1) - morrowair(x,h,2)), 1e7);

figure(1)
clf
loglog(e,alpha_air1,'r')
hold on
loglog(e,alpha_morrow,'b')
legend('air1', 'morrowair', 'Location', 'Northwest')
grid on 
title('Effective ionization coefficient \alpha_{eff}')
xlim([1e6 1e9])
xlabel('E (V/m)')
ylabel('\alpha_{eff} (m^{-1})')
set(gca, 'Fontsize', 15)

figure(2)
clf
loglog(e,nu_air1,'r')
hold on
loglog(e,nu_morrow,'b')
legend('air1', 'morrowair', 'Location', 'Northwest')
grid on 
title('Effective ionization frequency \nu_{eff}')
xlim([1e6 1e9])
set(gca, 'Fontsize', 15)





%%

%  clear all; clc;
% 
% e=10.^(5:0.00001:10);
% nui_air1=air1(e,0,10);
% nui_morrow=morrowair(e,0,1);
% nua_air1=air1(e,0,2);
% nua_morrow=morrowair(e,0,2);
% mue_air1 = air1(e,0,11);
% mue_morrow = morrowair(e,0,4);
% ve_air1 = -mue_air1.*e;
% ve_morrow = -mue_morrow.*e;
% alpha_i_air1 = nui_air1./abs(ve_air1);
% alpha_a_air1 = nua_air1./abs(ve_air1);
% 
% alpha_i_morrow = nui_morrow./abs(ve_morrow);
% alpha_a_morrow = nua_morrow./abs(ve_morrow);
% 
% figure(3)
% clf
% alpha_air1 = alpha_i_air1-alpha_a_air1;
% alpha_morrow = alpha_i_morrow-alpha_a_morrow;
% loglog(e,alpha_air1,'k')
% hold on
% loglog(e,alpha_morrow,'k')
% %plot(Ek_air1,alpha_i_air1(idx),'*k')
% %L=legend('\alpha_i - \alpha_{att} (air1)', '\alpha_i - \alpha_{att} (morrow)','Location', 'NW');
% xlabel('Electric field (V/m)')
% ylabel('\alpha_{ion} - \alpha_{att} (m^{-1})')
% %set(L,'FontSize',10);
% xlim([1e6 1e9])
% set(gca,'FontSize',15);
% grid on
% 

%%


R = 0.25e-3;  % [m]

eps0 = 8.85e-12; % vacuum permittivity
z = (0:0.01:0.5) * 1e-3;  % [m]


E0 =  [ 2      4      6      8    ] * 1e5;  % [v/m]
Qdr = [ 97.7,  91.0,  84.0,  76.7  ; ...
       238.0, 216.1, 193.6, 170.2 ] * 1e-12; % [C]
   
E = cell(2,4);

Alpha_i_air1 = cell(2,4);
Alpha_a_air1 = cell(2,4);

Alpha_i_morrow = cell(2,4);
Alpha_a_morrow = cell(2,4);



for i = 1:length(E0)   
    E{1,i} = E0(i) * ( 1+2*R^3./(R+z).^3 ) + Qdr(1,i) ./ ( 4*pi*eps0.*(R+z).^2 );  
    E{2,i} = E0(i) * ( 1+2*R^3./(R+z).^3 ) + Qdr(2,i) ./ ( 4*pi*eps0.*(R+z).^2 );  
end

for i = 1 : length(E0)
    [Alpha_i_air1{1,i}, Alpha_a_air1{1,i}] = CalculateAlpha(E{1,i}, h, 'air1');
    [Alpha_i_air1{2,i}, Alpha_a_air1{2,i}] = CalculateAlpha(E{2,i}, h, 'air1');
    
    [Alpha_i_morrow{1,i}, Alpha_a_morrow{1,i}] = CalculateAlpha(E{1,i}, h, 'morrowair');
    [Alpha_i_morrow{2,i}, Alpha_a_morrow{2,i}] = CalculateAlpha(E{2,i}, h, 'morrowair');
end




%% data point 1_1
figure(101)
clf
plot(z*1e3,E{1,1})
xlabel('z (mm)')
ylabel('E (V/m)')
set(gca, 'Fontsize', 15)

figure(102)
clf
loglog(E{1,1}, Alpha_i_air1{1,1} - Alpha_a_air1{1,1},'r')
hold on 
loglog(E{1,1}, Alpha_i_morrow{1,1} - Alpha_a_morrow{1,1},'b-.')
legend('air1', 'morrowair')
xlabel('E (V/m)')
ylabel('\alpha_{eff}')
title('point 1')
set(gca, 'Fontsize', 15)
saveas(gcf,'Point1_E_Alpha.png')

figure(103)
clf
semilogy(z*1e3, Alpha_i_air1{1,1} - Alpha_a_air1{1,1},'r')
hold on
semilogy(z*1e3, Alpha_i_morrow{1,1} - Alpha_a_morrow{1,1},'b-.')
legend('air1', 'morrowair')
xlabel('z (mm)')
ylabel('\alpha_{eff}')
title('point 1')
set(gca, 'Fontsize', 15)
saveas(gcf,'Point1_z_Alpha.png')


%% data point 1_2
figure(201)
clf
plot(z*1e3,E{1,2})
xlabel('z (mm)')
ylabel('E (V/m)')
set(gca, 'Fontsize', 15)

figure(202)
clf
loglog(E{1,2}, Alpha_i_air1{1,2} - Alpha_a_air1{1,2},'r')
hold on 
loglog(E{1,2}, Alpha_i_morrow{1,2} - Alpha_a_morrow{1,2},'b-.')
legend('air1', 'morrowair')
xlabel('E (V/m)')
ylabel('\alpha_{eff}')
title('point 2')
set(gca, 'Fontsize', 15)
saveas(gcf,'Point2_E_Alpha.png')

figure(203)
clf
semilogy(z*1e3, Alpha_i_air1{1,2} - Alpha_a_air1{1,2},'r')
hold on
semilogy(z*1e3, Alpha_i_morrow{1,2} - Alpha_a_morrow{1,2},'b-.')
legend('air1', 'morrowair')
xlabel('z (mm)')
ylabel('\alpha_{eff}')
title('point 2')
set(gca, 'Fontsize', 15)
saveas(gcf,'Point2_z_Alpha.png')

%% data point 1_3
figure(201)
clf
plot(z*1e3,E{1,3})
xlabel('z (mm)')
ylabel('E (V/m)')
set(gca, 'Fontsize', 15)

figure(202)
clf
loglog(E{1,3}, Alpha_i_air1{1,3} - Alpha_a_air1{1,3},'r')
hold on 
loglog(E{1,3}, Alpha_i_morrow{1,3} - Alpha_a_morrow{1,3},'b-.')
legend('air1', 'morrowair')
xlabel('E (V/m)')
ylabel('\alpha_{eff}')
title('point 3')
set(gca, 'Fontsize', 15)
saveas(gcf,'Point3_E_Alpha.png')

figure(203)
clf
semilogy(z*1e3, Alpha_i_air1{1,3} - Alpha_a_air1{1,3},'r')
hold on
semilogy(z*1e3, Alpha_i_morrow{1,3} - Alpha_a_morrow{1,3},'b-.')
legend('air1', 'morrowair')
xlabel('z (mm)')
ylabel('\alpha_{eff}')
title('point 3')
set(gca, 'Fontsize', 15)
saveas(gcf,'Point3_z_Alpha.png')


%% data point 1_4
figure(201)
clf
plot(z*1e3,E{1,4})
xlabel('z (mm)')
ylabel('E (V/m)')
set(gca, 'Fontsize', 15)

figure(202)
clf
semilogx(E{1,4}, Alpha_i_air1{1,4} - Alpha_a_air1{1,4},'r')
hold on 
semilogx(E{1,4}, Alpha_i_morrow{1,4} - Alpha_a_morrow{1,4},'b-.')
legend('air1', 'morrowair')
xlabel('E (V/m)')
ylabel('\alpha_{eff}')
title('point 4')
set(gca, 'Fontsize', 15)
saveas(gcf,'Point4_E_Alpha.png')

figure(203)
clf
plot(z*1e3, Alpha_i_air1{1,4} - Alpha_a_air1{1,4},'r')
hold on
plot(z*1e3, Alpha_i_morrow{1,4} - Alpha_a_morrow{1,4},'b-.')
legend('air1', 'morrowair')
xlabel('z (mm)')
ylabel('\alpha_{eff}')
title('point 4')
set(gca, 'Fontsize', 15)
saveas(gcf,'Point4_z_Alpha.png')


%% data point 2_1
figure(701)
clf
plot(z*1e3,E{2,1})
xlabel('z (mm)')
ylabel('E (V/m)')
set(gca, 'Fontsize', 15)

figure(702)
clf
plot(E{2,1}, Alpha_i_air1{2,1} - Alpha_a_air1{2,1},'r')
hold on 
plot(E{2,1}, Alpha_i_morrow{2,1} - Alpha_a_morrow{2,1},'b-.')
legend('air1', 'morrowair')
xlabel('E (V/m)')
ylabel('\alpha_{eff}')
set(gca, 'Fontsize', 15)

figure(703)
clf
plot(z*1e3, Alpha_i_air1{2,1} - Alpha_a_air1{2,1},'r')
hold on
plot(z*1e3, Alpha_i_morrow{2,1} - Alpha_a_morrow{2,1},'b-.')
legend('air1', 'morrowair')
xlabel('z (mm)')
ylabel('\alpha_{eff}')
set(gca, 'Fontsize', 15)


%% data point 2_2
figure(801)
clf
plot(z*1e3,E{2,2})
xlabel('z (mm)')
ylabel('E (V/m)')
set(gca, 'Fontsize', 15)

figure(802)
clf
plot(E{2,2}, Alpha_i_air1{2,2} - Alpha_a_air1{2,2},'r')
hold on 
plot(E{2,2}, Alpha_i_morrow{2,2} - Alpha_a_morrow{2,2},'b-.')
legend('air1', 'morrowair')
xlabel('E (V/m)')
ylabel('\alpha_{eff}')
set(gca, 'Fontsize', 15)

figure(803)
clf
plot(z*1e3, Alpha_i_air1{2,2} - Alpha_a_air1{2,2},'r')
hold on
plot(z*1e3, Alpha_i_morrow{2,2} - Alpha_a_morrow{2,2},'b-.')
legend('air1', 'morrowair')
xlabel('z (mm)')
ylabel('\alpha_{eff}')
set(gca, 'Fontsize', 15)


%% data point 2_3
figure(701)
clf
plot(z*1e3,E{2,3})
xlabel('z (mm)')
ylabel('E (V/m)')
set(gca, 'Fontsize', 15)

figure(702)
clf
plot(E{2,3}, Alpha_i_air1{2,3} - Alpha_a_air1{2,3},'r')
hold on 
plot(E{2,3}, Alpha_i_morrow{2,3} - Alpha_a_morrow{2,3},'b-.')
legend('air1', 'morrowair')
xlabel('E (V/m)')
ylabel('\alpha_{eff}')
set(gca, 'Fontsize', 15)

figure(703)
clf
plot(z*1e3, Alpha_i_air1{2,3} - Alpha_a_air1{2,3},'r')
hold on
plot(z*1e3, Alpha_i_morrow{2,3} - Alpha_a_morrow{2,3},'b-.')
legend('air1', 'morrowair')
xlabel('z (mm)')
ylabel('\alpha_{eff}')
set(gca, 'Fontsize', 15)


%% data point 2_4
figure(801)
clf
plot(z*1e3,E{2,4})
xlabel('z (mm)')
ylabel('E (V/m)')
set(gca, 'Fontsize', 15)

figure(802)
clf
plot(E{2,4}, Alpha_i_air1{2,4} - Alpha_a_air1{2,4},'r')
hold on 
plot(E{2,4}, Alpha_i_morrow{2,4} - Alpha_a_morrow{2,4},'b-.')
legend('air1', 'morrowair')
xlabel('E (V/m)')
ylabel('\alpha_{eff}')
set(gca, 'Fontsize', 15)

figure(803)
clf
plot(z*1e3, Alpha_i_air1{2,4} - Alpha_a_air1{2,4},'r')
hold on
plot(z*1e3, Alpha_i_morrow{2,4} - Alpha_a_morrow{2,4},'b-.')
legend('air1', 'morrowair')
xlabel('z (mm)')
ylabel('\alpha_{eff}')
set(gca, 'Fontsize', 15)


% %% data point 3_1
% figure(701)
% clf
% plot(z*1e3,E{3,1})
% xlabel('z (mm)')
% ylabel('E (V/m)')
% set(gca, 'Fontsize', 15)
% 
% figure(702)
% clf
% plot(E{3,1}, Alpha_i_air1{3,1} - Alpha_a_air1{3,1},'r')
% hold on 
% plot(E{3,1}, Alpha_i_morrow{3,1} - Alpha_a_morrow{3,1},'b-.')
% legend('air1', 'morrowair')
% xlabel('E (V/m)')
% ylabel('\alpha_{eff}')
% set(gca, 'Fontsize', 15)
% 
% figure(703)
% clf
% plot(z*1e3, Alpha_i_air1{3,1} - Alpha_a_air1{3,1},'r')
% hold on
% plot(z*1e3, Alpha_i_morrow{3,1} - Alpha_a_morrow{3,1},'b-.')
% legend('air1', 'morrowair')
% xlabel('z (mm)')
% ylabel('\alpha_{eff}')
% set(gca, 'Fontsize', 15)
% 
% 
% %% data point 2_2
% figure(801)
% clf
% plot(z*1e3,E{3,2})
% xlabel('z (mm)')
% ylabel('E (V/m)')
% set(gca, 'Fontsize', 15)
% 
% figure(802)
% clf
% plot(E{3,2}, Alpha_i_air1{3,2} - Alpha_a_air1{3,2},'r')
% hold on 
% plot(E{3,2}, Alpha_i_morrow{3,2} - Alpha_a_morrow{3,2},'b-.')
% legend('air1', 'morrowair')
% xlabel('E (V/m)')
% ylabel('\alpha_{eff}')
% set(gca, 'Fontsize', 15)
% 
% figure(803)
% clf
% plot(z*1e3, Alpha_i_air1{2,2} - Alpha_a_air1{2,2},'r')
% hold on
% plot(z*1e3, Alpha_i_morrow{2,2} - Alpha_a_morrow{2,2},'b-.')
% legend('air1', 'morrowair')
% xlabel('z (mm)')
% ylabel('\alpha_{eff}')
% set(gca, 'Fontsize', 15)
% 
% 
% %% data point 2_3
% figure(701)
% clf
% plot(z*1e3,E{2,3})
% xlabel('z (mm)')
% ylabel('E (V/m)')
% set(gca, 'Fontsize', 15)
% 
% figure(702)
% clf
% plot(E{2,3}, Alpha_i_air1{2,3} - Alpha_a_air1{2,3},'r')
% hold on 
% plot(E{2,3}, Alpha_i_morrow{2,3} - Alpha_a_morrow{2,3},'b-.')
% legend('air1', 'morrowair')
% xlabel('E (V/m)')
% ylabel('\alpha_{eff}')
% set(gca, 'Fontsize', 15)
% 
% figure(703)
% clf
% plot(z*1e3, Alpha_i_air1{2,3} - Alpha_a_air1{2,3},'r')
% hold on
% plot(z*1e3, Alpha_i_morrow{2,3} - Alpha_a_morrow{2,3},'b-.')
% legend('air1', 'morrowair')
% xlabel('z (mm)')
% ylabel('\alpha_{eff}')
% set(gca, 'Fontsize', 15)
% 
% 
% %% data point 2_4
% figure(801)
% clf
% plot(z*1e3,E{2,4})
% xlabel('z (mm)')
% ylabel('E (V/m)')
% set(gca, 'Fontsize', 15)
% 
% figure(802)
% clf
% plot(E{2,4}, Alpha_i_air1{2,4} - Alpha_a_air1{2,4},'r')
% hold on 
% plot(E{2,4}, Alpha_i_morrow{2,4} - Alpha_a_morrow{2,4},'b-.')
% legend('air1', 'morrowair')
% xlabel('E (V/m)')
% ylabel('\alpha_{eff}')
% set(gca, 'Fontsize', 15)
% 
% figure(803)
% clf
% plot(z*1e3, Alpha_i_air1{2,4} - Alpha_a_air1{2,4},'r')
% hold on
% plot(z*1e3, Alpha_i_morrow{2,4} - Alpha_a_morrow{2,4},'b-.')
% legend('air1', 'morrowair')
% xlabel('z (mm)')
% ylabel('\alpha_{eff}')
% set(gca, 'Fontsize', 15)








   

















