clear; close all; clc;


R = [0.75e-3, 0.5e-3, 0.25e-3]; % m
E0 = [2e5, 4e5, 6e5, 8e5]; % v/m

QdrMat = zeros(length(E0), length(R));
Meek0Mat = zeros(length(E0), length(R));

QdrMatmorrow = zeros(length(E0), length(R));
Meek0Matmorrow = zeros(length(E0), length(R));

h = fzero( @(h) satm(h)-0.4, 8); % find the altitude corresponding to 0.4 atm



for i = 1:length(E0)
    for j = 1:length(R)
         QdrMat(i,j) = fzero(@(x) oneSphere(x, E0(i), R(j), h, 'air1'), 1e-12);  % x is corresponding to Qdr
%         QdrMat(i,j) = fzero(@(x) oneSphere(x, E0(i), R(j), h, 'air1'), 3e-10);
         Meek0Mat(i,j) = oneSphere(QdrMat(i,j) , E0(i), R(j),  h, 'air1'); % store values, check whether Meek=0 is satisfied
         QdrMatmorrow(i,j) = fzero(@(x) oneSphere(x, E0(i), R(j), h, 'morrowair'), 3e-10);
         Meek0Matmorrow(i,j) = oneSphere(QdrMatmorrow(i,j) , E0(i), R(j),  h, 'morrowair'); % store values, check whether Meek=0 is satisfied 
    end
end



% Qdr[pC] values read from Babitch
q1B = [400, 355, 315, 270]; % R = 0.75mm
q2B = [220, 200, 178, 157]; % R = 0.5mm
q3B = [80, 75, 69, 63]; % R = 0.25mm

% Qdr [pC] values from looking tables
q1 = 395:-43:266;
q2 = 213:-20:153;
q3 = 76:-5:61;



figure(1)
clf
hold on
plot(E0, QdrMat(:,1)*1e12,'o--') % convert from C to pC
plot(E0, QdrMat(:,2)*1e12,'*--')
plot(E0, QdrMat(:,3)*1e12,'s--')

plot(E0, QdrMatmorrow(:,1)*1e12,'o-.') % convert from C to pC
plot(E0, QdrMatmorrow(:,2)*1e12,'*-.')
plot(E0, QdrMatmorrow(:,3)*1e12,'s-.')

plot(E0,q1B,'o-')
plot(E0,q2B,'*-')
plot(E0,q3B,'s-')

% figure(1)
% clf
% hold on
% leg1 = plot(E0*1e-5, QdrMat(:,1)*1e12,'k*--'); % convert from v/m to kv/cm; C to pC
% plot(E0*1e-5, QdrMat(:,2)*1e12,'k*--')
% plot(E0*1e-5, QdrMat(:,3)*1e12,'k*--')
% 
% leg2 = plot(E0*1e-5, QdrMatmorrow(:,1)*1e12,'kx-.'); % convert from C to pC
% plot(E0*1e-5, QdrMatmorrow(:,2)*1e12,'kx-.')
% plot(E0*1e-5, QdrMatmorrow(:,3)*1e12,'kx-.')
% 
% leg3 = plot(E0*1e-5,q1B,'ko-');
% plot(E0*1e-5,q2B,'ko-')
% plot(E0*1e-5,q3B,'ko-')

%legend([leg1, leg2, leg3],  {'Meek=20 (air1)', 'Meek=20 (morrowair)',  'Babich et al., 2016'}); 


L=legend('R=0.75mm Meek (air1)','R=0.5mm Meek (air1)','R=0.25mm Meek (air1)','R=0.75mm Meek (morrow)','R=0.5mm Meek (morrow)','R=0.25mm Meek (morrow)', 'R=0.75mm Babitch','R=0.5mm Babitch','R=0.25mm Babitch', 'FontName', 'helvetica');
% xlim([1,8.5])
% ylim([25, 450])
xlabel('Applied Electric Field (kV/cm)')
ylabel('Q_{dr} (pC)')
% grid on
%title('Q_{dr} as a function of aplied electric field (fzero approach)')

% set(L,'FontName','times','FontSize',12);
set(gca,'fontsize', 15)
%set(gcf,'PaperPositionMode','auto');


return
% Qdr [pC] values from looking tables
q1 = 395:-43:266;
q2 = 213:-20:153;
q3 = 76:-5:61;


figure(2)
hold on
plot(E0, q1,'o-') % convert from C to pC
plot(E0, q2,'*-')
plot(E0, q3,'s-')

plot(E0,q1B,'o--')
plot(E0,q2B,'*--')
plot(E0,q3B,'x--')


L=legend('R=0.75mm Meek','R=0.5mm Meek','R=0.25mm Meek','R=0.75mm Babitch','R=0.5mm Babitch','R=0.25mm Babitch');

xlabel('Applied Electric Field (kV/cm)')
ylabel('Q_{dr} (pC)')
grid on
title('Q_{dr} as a function of aplied electric field (lookup table approach)')

set(L,'FontName','times','FontSize',12);
set(gca,'FontName','times','FontSize',12);
set(gcf,'PaperPositionMode','auto');
