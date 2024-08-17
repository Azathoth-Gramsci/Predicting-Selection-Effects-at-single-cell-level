function n = state2n2dPlus(z,xmax,s)
    
    %Computes n from the state
    if z(3) == 1
        n = z(2)*(xmax+1) + z(1) + 1;
    else
        n = z(2)*(xmax+1) + z(1) + 1 + floor(s/2);
    end

end