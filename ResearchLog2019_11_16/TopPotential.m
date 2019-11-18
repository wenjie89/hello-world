function potential = TopPotential(z,r,Z,R,S)
% This function cauculates the potential (only one point) at the top
% boundary due to the source distribution.
% z: z coordinates of the point for the potential on the right boundary
% r: r coordinates of the point for the potential on the right boundary
% Z: Z coordinates of the whole computational domain after meshgrid
% R: R coordinates of the whole computational domain after meshgrid


[idxx,idxy] = find(S~=0);


eps0 =  8.85e-12;
phi = 0;
% phi_test = 0;

for i = 1:length(idxx)
    xgoal = idxx(i); ygoal = idxy(i);  % get the x y index coor. to source at i idx
    
    r_source = R(xgoal,ygoal); z_source =  Z(xgoal,ygoal);  % r and z of the source that leads to the potential
    
    phi = phi + 1/(4*pi*eps0)*S( xgoal,ygoal) / sqrt((r_source-r)^2 + (z_source-z)^2);
    
   % phi_test = phi_test + 1/(4*pi*eps0)*S(idxx(i),idxy(i)) / sqrt((R(idxx(i),idxy(i))-roro)^2 + (Z(idxx(i),idxy(i))-zz)^2);
    
end

potential = phi;

end


