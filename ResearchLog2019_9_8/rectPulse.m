function p = rectPulse(x)

    p = ones(1,length(x));
    
    idx1 = find(x<0);
    idx2 = find(x>1);
    
    p(idx1) = 0;
    p(idx2) = 0;
    
end
