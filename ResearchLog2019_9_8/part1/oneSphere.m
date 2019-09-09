function Meek0 = oneSphere(Qdr, E0, R, h, option)

%function [alpha_i, alpha_a, alpha, Meek0] = oneSphere(Qdr, E0, R)
    e0 = 8.85e-12; % vacuum permittivity

    dz = 1e-6;

    z = 0:dz:0.01;

    % total E field
    E = E0*(1+2*R^3./(R+z).^3)+Qdr./(4*pi*e0*(R+z).^2);
    
    % nu
    if strcmp(option, 'air1')
        nui = air1(E,h,10);
        nua = air1(E,h,2);
        mue = air1(E,h,11);
    elseif strcmp(option, 'morrowair')
        nui = morrowair(E,h,1);
        nua = morrowair(E,h,2);
        mue = morrowair(E,h,4);
    else 
        error('Wrong option!')
    end
    
   
%         figure
%         plot(z,nui-nua)
%         find(abs(nui-nua)<0.1,1)
%     
%        length( find(abs(nui-nua)<10) )
        
    ve = -mue.*E;
    
    % alpha
    alpha_i = nui./abs(ve);
    alpha_a = nua./abs(ve);
    alpha = alpha_i-alpha_a;
    
%      if find(alpha<0)==0
%         error('Increase the length of simulation to find zin !')
%      end
%      if strcmp(option, 'air1')
%         zin = fzero( @(x) air1(E0*(1+2*R^3./(R+x).^3)+Qdr./(4*pi*e0*(R+x).^2),h,10) - ...
%                      air1(E0*(1+2*R^3./(R+x).^3)+Qdr./(4*pi*e0*(R+x).^2),h,2), z(find(alpha<0,1)) );
%      elseif strcmp(option, 'morrowair')
%          zin = fzero( @(x) morrowair(E0*(1+2*R^3./(R+x).^3)+Qdr./(4*pi*e0*(R+x).^2),h,1) - ...
%                      air1(E0*(1+2*R^3./(R+x).^3)+Qdr./(4*pi*e0*(R+x).^2),h,2), z(find(alpha<0,1)) );
%      end
                 
                 
%     alphaPos = alpha(alpha>=0);
%     n = numel(alphaPos);
%     fprintf('n=%d\n',n);
   
    alpha(alpha<0) = 0;
%     m = numel(alpha);
%     fprintf('m=%d\n',m);
    
    % Meek 
    Meek0 = sum(alpha) * dz - 20;
    
    
end



    
