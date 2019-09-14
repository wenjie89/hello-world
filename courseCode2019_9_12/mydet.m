function detA = mydet(A)

    [~, U, ~, sig] = lutx(A);  % First, perform LU decomposition on A. 
                               % L is not needed to be computed since we already
                               % know that its determinant is one. P is not
                               % needed as well, since the determinant of P
                               % is only related to the number of permutations.

    detU = prod(diag(U));      % the determinant of U

    detP = sig;                % the determinant of P 
    
    detA = detU/detP;          % det(L)det(U) = det(P)det(A), det(L) = 1, 
                               % therefore det(A) = det(U)/det(P)  

end











