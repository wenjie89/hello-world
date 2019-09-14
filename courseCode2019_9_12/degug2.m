%% original script

close all; clear; clc


A = [1.5 0.8 5.1; 3 4 5; -3 -2 -6];
Aori = A;



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
             kshow1 = k
             Ainterm1 = A
             p([k m]) = p([m k]);
             numPiv = numPiv+1; % increase the number by 1 while pivoting
        end
     
        % Compute multipliers
        Ainterm12 = A;
        i = k+1:n
        kshow2 = k
        Aik = A(i,k)  % before processing
        Akk = A(k,k)  % before processing
        A(i,k) = A(i,k)/A(k,k);  % processing
        Aikfin = A(i,k)          % after processing
       
        Ainterm2 = A             % after processing
        % Update the remainder of the matrix
        j = k+1:n
        i
        A(i,j) = A(i,j) - A(i,k)*A(k,j);
   
        Ainterm3 = A
   end
end
% Separate result
L = tril(A,-1) + eye(n,n);
U = triu(A);


    