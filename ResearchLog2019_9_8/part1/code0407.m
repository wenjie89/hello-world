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


figure(1)
clf
loglog(e,alpha_i_air1,'k')
hold on
loglog(e,alpha_i_morrow,'k')
loglog(e,alpha_a_air1,'k--')
loglog(e,alpha_a_morrow,'k--')
% idx = find(alpha_i_air1-alpha_a_air1>0,1);
% Ek0 = e(idx); % initial guess of Ek
% Ek_air1 = fzero(@(e) air1(e,0,10) - air1(e,0,2), Ek0); % find the exact Ek 
% Ek_morrow = fzero(@(e) morrowair(e,0,1) - morrowair(e,0,2), Ek0); % find the exact Ek 
% plot(Ek_air1,alpha_i_air1(idx),'*k')
%L=legend('\alpha_i (air1)','\alpha_i (morrow)','\alpha_{att} (air1)', '\alpha_{att} (morrow)', 'Location', 'NW');
xlabel('Electric field (V/m)')
ylabel('Ionization and Attachment coefficient (m^{-1})')
xlim([1e6 1e9])
%set(L,'FontSize',10);
set(gca,'FontSize',15)
%set(gcf,'PaperPositionMode','auto');


figure(2)
clf
loglog(e,ve_air1,'k')
hold on
loglog(e,ve_morrow,'k')
xlabel('Electric field (V/m)')
ylabel('v_e (m/s)')
%L=legend('v_e (air1)','v_e (morrow) ');
%set(L,'FontSize',10);
xlim([1e6 1e9])
set(gca,'FontSize',15)

%set(gcf,'PaperPositionMode','auto');

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

%set(gcf,'PaperPositionMode','auto');



figure(4)
clf
loglog(e, nui_air1, 'k-', e, nua_air1, 'k--')
hold on 
loglog(e, nui_morrow, 'k-', e, nua_morrow, 'k--')
%plot(Ek_air1,nui_air1(idx),'*k')
%axis([1e5 1e7 1e7 1e11])
%L=legend('\nu_i (air1)','\nu_{att} (air1)','\nu_i (morrow)','\nu_{att} (morrow)','location','NW');
xlabel('Electric field (V/m)')
ylabel('Ionization and Attachment frequency (s^{-1})')
xlim([1e6 1e9])
%set(L,'FontSize',10);
set(gca,'FontSize',15);
%set(gcf,'PaperPositionMode','auto');




