function p = rectPulse(x)

    p = ones(1,length(x));
    delx = x(2) - x(1);
    % idx1 = find(x<0);
    idx1 = find(0-x > 0.5*delx);
    % idx2 = find(x>1);
    idx2 = find(x-1 > 0.5*delx);
    
    p(idx1) = 0;
    p(idx2) = 0;
    
end
