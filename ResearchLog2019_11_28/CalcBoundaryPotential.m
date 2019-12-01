function potential = CalcBoundaryPotential(z,r,Z,R,S,V)
% This function cauculates the potential (only one point) at the boundary 
% due to the source distribution.
% z: z coordinates of the point for the potential on the boundary
% r: r coordinates of the point for the potential on the boundary
% Z: Z coordinates of the whole computational domain after meshgrid
% R: R coordinates of the whole computational domain after meshgrid
% V: volume discretization


[idxx,idxy] = find(S~=0);


eps0 =  8.85e-12;
phi = 0;
% phi_test = 0;

for i = 1:length(idxx)
    xgoal = idxx(i); ygoal = idxy(i);  % get the x y index coor. to source at i idx
    
    r_source = R(xgoal,ygoal); z_source =  Z(xgoal,ygoal);  % r and z of the source that leads to the potential
    
   % phi = phi + 1/(4*pi*eps0)*S( xgoal,ygoal)*V(xgoal,ygoal) / sqrt((r_source-r)^2 + (z_source-z)^2);

    phi = phi + 1/(4*pi*eps0)*S( xgoal,ygoal)*V(xgoal,ygoal) / sqrt((r_source+r)^2 + (z_source-z)^2) * ...
        4 * Kappa(sqrt(4*r*r_source/((r_source+r)^2 + (z_source-z)^2))) / (2*pi);
    
    
    
   % phi_test = phi_test + 1/(4*pi*eps0)*S(idxx(i),idxy(i)) / sqrt((R(idxx(i),idxy(i))-roro)^2 + (Z(idxx(i),idxy(i))-zz)^2);
    
end

potential = phi;

end


function kappa = Kappa(k)

t=1-k*k;

kappa=(((0.032024666*t+0.054544409)*t + 0.097932891)*t + 1.3862944)...
- (((0.010944912*t +0.060118519)*t+0.12475074)*t+0.5)*log(t);

end

