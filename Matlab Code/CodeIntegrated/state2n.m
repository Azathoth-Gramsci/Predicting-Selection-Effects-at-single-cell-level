function n = state2n(z,xmax,ymax)
    global s
    
    %Computes n from the state
    n = z(2)*(xmax+1) + z(1) + 1;

    %All indices outside the truncation are part of J' and 
    %are put at the end of their corresponding row.
    if z(1)>xmax || z(2)>ymax 
        n = s;
    end
end