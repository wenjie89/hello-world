close all; clear all; clc;

e=10.^(5:0.00001:10);
nui_air1=air1(e,0,10);
nui_morrow=morrowair(e,0,1);
nua_air1=air1(e,0,2);
nua_morrow=morrowair(e,0,2);
mue_air1 = air1(e,0,11);
mue_morrow = morrowair(e,0,4);
ve_air1 = -mue_air1.*e;
ve_morrow = -mue_morrow.*e;
alpha_i_air1 = nui_air1./abs(ve_air1);
alpha_a_air1 = nua_air1./abs(ve_air1);

alpha_i_morrow = nui_morrow./abs(ve_morrow);
alpha_a_morrow = nua_morrow./abs(ve_morrow);


figure(3)
clf
alpha_air1 = alpha_i_air1-alpha_a_air1;
alpha_morrow = alpha_i_morrow-alpha_a_morrow;
loglog(e,alpha_air1,'k')
hold on
loglog(e,alpha_morrow,'k')
%plot(Ek_air1,alpha_i_air1(idx),'*k')
%L=legend('\alpha_i - \alpha_{att} (air1)', '\alpha_i - \alpha_{att} (morrow)','Location', 'NW');
xlabel('Electric field (V/m)')
ylabel('\alpha_{ion} - \alpha_{att} (m^{-1})')
%set(L,'FontSize',10);
xlim([1e6 1e9])
set(gca,'FontSize',15);