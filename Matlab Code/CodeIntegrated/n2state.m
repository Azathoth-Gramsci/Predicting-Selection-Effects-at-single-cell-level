function z = n2state(n,xmax)
%Computes the state from n
x = mod(n-1,xmax+1);  %can be zero
y = (n-1-x)/(xmax+1); %can be zero
z = [x,y]';
end