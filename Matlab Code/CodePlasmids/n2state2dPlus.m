function z = n2state2dPlus(n,xmax,s)
    %Computes the state from n
    if n < s/2
        x = mod(n-1,xmax+1);
        y = (n-1-x)/(xmax+1);
        z = [x,y,1]';
    else
        x = mod(n-floor(s/2)-1,xmax+1);
        y = (n-floor(s/2)-1-x)/(xmax+1);
        z = [x,y,2]';
    end
end