function A = addprop(z,n,ni,A,mu,xmax,s)
    %Adds propensities to the reaction matrix
    global sysP
    
    if mu == 1 %bursty production
       burst = sysP.rcp(3);
       pAdd = 0;
       x = z(1);
       for xnew = x+1:xmax-1
           znew = z;
           znew(1) = xnew;
           nnew = state2n(znew,xmax,2);
           A(n,nnew) = A(n,nnew) + a1(z)*(1-1/burst)^(xnew-x-1)*1/burst; %redistribute the rate a according to a geometric distribution;
           pAdd = pAdd + (1-1/burst)^(xnew-x-1)*1/burst; %save how much prob has been already redistributed
       end
       xnew = xmax;
       znew = z;
       znew(1) = xnew;
       nnew = state2n(znew,xmax,2);
       A(n,nnew) = a1(z)*(1-pAdd); %assign all mass that would go out of the truncation are to the last truncated state
    elseif mu == 2 %protein dilution
       A(n,ni) = A(n,ni) + a2(z);
    elseif mu == 3 %differentiation
       A(n,ni) = A(n,ni) + a3(z);
    end
end

%PROPENSITIES 

function out = a1(z)
global sysP
a = sysP.rcp(1);
out = a;
end

function out = a2(z)
global sysP
b = sysP.rcp(2);
out = b*z(1);
end


function out = a3(z)
global sysP
rH = sysP.rcp(4); 
Kh = sysP.rcp(5); 
nH = sysP.rcp(6);
uIn = sysP.rcp(7);
parasDiffFunction = [rH,Kh,nH];
out = uIn*DifferentiationFunction(z(1),parasDiffFunction)*(1-z(2));
end

