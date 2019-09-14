close all; clear; clc;

A = [1.5 0.8 5.1; 3 4 5; -3 -2 -6];


% [L, U, p, sig] = lutx(A);
% [Lmod, Umod, pmod, sigmod] = lutxMod(A);
% [Lmat, Umat, Pmat] = lu(A);







% function [L, U, p, sig] = lutx(A)
%LU Triangular factorization
%   [L,U,p] = lutx(A) produces a unit lower triangular
%   matrix L, an upper triangular matrix U, and a
%   permutation vector p, so that L*U = A(p,:).

numPiv = 0;  % initialize the number of pivots

[n,n] = size(A);
p = (1:n)';

for k = 1:n-1
   % Find largest element below diagonal in k-th column
   [r,m] = max(abs(A(k:n,k)));
   m = m+k-1;
   % Skip elimination if column is zero
   if (A(m,k)  ~= 0)
        % Swap pivot row
        if (m ~= k)
             A([k m],:) = A([m k],:);
             p([k m]) = p([m k]);
             numPiv = numPiv+1; % increase the number by 1 while pivoting
        end
        % Compute multipliers
        i = k+1:n;
        A(i,k) = A(i,k)/A(k,k);
        % Update the remainder of the matrix
        j = k+1:n;
        A(i,j) = A(i,j) - A(i,k)*A(k,j);
   end
end
% Separate result
L = tril(A,-1) + eye(n,n);
U = triu(A);

% calculate sig
if rem(numPiv,2) == 0      % even number of permutation
    sig = 1;
elseif rem(numPiv,2) == 1  % odd number of permutation
    sig = -1;
end
    

% end                   
            
        
  
    
    
%%
Anew = A;

A = [1.5 0.8 5.1; 3 4 5; -3 -2 -6];


[Lori, Uori, pori, sigori] = lutx(A);
[Lmod, Umod, pmod, sigmod] = lutxMod(A);
[Lmat, Umat, Pmat] = lu(A);


                   
[L, U, p, sig, B] = lutx(A);
[Lmod, Umod, pmod, sigmod, Bmod] = lutxMod(A);        
           
           
     
        
        
       
      
       
       
       
    
        
      


