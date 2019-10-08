function [f,it] = mySOR_cylinder(f,r,z,rou,w,epsilon)
% f --       n x m, electric potential 
% r --       1 x m, radius vector
% z --       1 x n, longitudinal vector
% rou --     n x m, source distribution
% w --       1 x 1, relaxation factor
% epsilon -- 1 x 1, accuracy factor



% EE 430 Spring 2019
%
% [f,it]=sorEE430S18(a,b,c,d,f,s,w,epsilon,nx,ny,dnorms) solves
% five-point difference equation for f(ix,iy) with constant coefficents
% a,b,c and d:
%
% a*f(ix,iy-1)+b*f(ix,iy+1)+c*f(ix-1,iy)+d*f(ix+1,iy)+f(ix,iy)=s(ix,iy)
%
% The solution is obtained in a domain consisting of nx+1 by ny+1 grid ponts.
% The boundary values of f (i.e., at ix=1, ix=nx+1, iy=1:ny+1; and iy=1, iy=ny+1, ix=1:nx+1)
% should be specified in matrix f before calling the function, 
% and values of f inside of the domain (ix=2:nx, iy=2:ny) 
% are obtained by successive overrelaxation (SOR) technique.
%
% w is a relaxation factor (1<w<2). The best value is 
% w=2/(1+sqrt(1-rosor^2)), where in the case of 
% the two-dimensional Poisson equation in the square (nx=ny) 
% rosor=cos(pi/nx).
%
% epsilon - desired accuracy of the solution.
%
% A good approximation for dnorms is
%
% dnorms=sum(sum(s.^2))
%
% Reference: 
%
% Hockney, R. W. and J. W. Eastwood, Computer simulation using particles,
% McGraw-Hill, 1981, pages 174-181


    % define parameters-------------------------------------------------------------------
    eps0 = 8.85e-12;
    delr = r(2)-r(1);
    delz = z(2)-z(1);

    u = -(2/delr^2+2/delz^2);                u_ = -(4/delr^2+2/delz^2); 
    a = (-1./(2*r*delr)+1/delr^2)/u;
    b = ( 1./(2*r*delr)+1/delr^2)/u;         b_ = 4/delr^2/u_;
    c = 1/delz^2/u;                          c_ = 1/delz^2/u_;
    d = c;                                   d_ = c_;
    s = -rou/eps0/u;                         s_ = -rou(:,1)/eps0/u_;          s(:,1) = s_;


    itMax = 30000;
    it = 0;
    norm_s = sum(sum(s.^2));

    norm_f_store = [];
    
    while (1)
        it = it+1;
        norm_r = 0;
        
        
        
        
        % main algorithm------------------------------------------------------------------

        % main part 
        for i = 2:size(f,1)-1
            %fprintf('i=%d\n',i);
            for j = 2:size(f,2)-1
                
                %fprintf('j=%d\n',j);
                r = a(j)*f(i,j-1)+b(j)*f(i,j+1)+c*f(i-1,j)+d*f(i+1,j)+f(i,j)-s(i,j);
                f(i,j) = f(i,j) - r * w;
                norm_r = norm_r + r^2;
            end
        end

        % boundary
        for i = 2:size(f,1)-2
            r_ = b_*f(i,2)+c_*f(i-1,1)+d_*f(i+1,1)+f(i,1)-s(i,1);
            f(i,1) = f(i,1) - r_ * w;
            norm_r = norm_r + r_^2;
        end
        
        
        
        % stopping condition--------------------------------------------------------------
        
        % condition (1)



       

        % condition (2)
        if  norm_r/norm_s < epsilon    % precision requirement is satisfied
            break;
        end


        % condition (3)
        if it > itMax 
             fprintf('The solver does not converge after %d iterations.\n', itMax);
            break;
           
        end

    end




end 
