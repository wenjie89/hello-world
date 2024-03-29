function [f,it,dnorm]=sorEE430S19(a,b,c,d,f,s,w,epsilon,nz,nr,dnorms)
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


it=0; %number of iterations

while(1)

it=it+1;
dnorm=0;
for iz=2:nz
    for jr=2:nr
        residual=a(jr)*f(iz,jr-1)+b(jr)*f(iz,jr+1)+...
        c*f(iz-1,jr)+d*f(iz+1,jr)+f(iz,jr)-s(iz,jr);
        dnorm=dnorm+residual^2;
        f(iz,jr)=f(iz,jr)-w*residual;
    end%iy
end%ix

if dnorm/dnorms < epsilon 
	break 
end%if

if it> 10000
    warning(' Solution has not been reached in 10000 iterations.')
    break
end%if

end %while
