function A = addpropNoDiff(z,n,ni,A,mu,xmax,ymax)
    %Adds propensities to the reaction matrix
    global sysP
    
    if mu == 1 %bursty production
           burst = sysP.rcp(3);
           pAdd = 0;
           x = z(1);
           for xnew = x+1:xmax-1
               znew = z;
               znew(1) = xnew;
               nnew = state2n(znew,xmax,ymax);
               A(n,nnew) = A(n,nnew) + a1(z)*(1-1/burst)^(xnew-x-1)*1/burst; %redistribute the rate a according to a geometric distribution;
               pAdd = pAdd + (1-1/burst)^(xnew-x-1)*1/burst; %save how much prob has been already redistributed
           end
           xnew = xmax;
           znew = z;
           znew(1) = xnew;
           nnew = state2n(znew,xmax,ymax);
           A(n,nnew) = a1(z)*(1-pAdd); %assign all mass that would go out of the truncation are to the last truncated state
    elseif mu == 2 %protein dilution
           A(n,ni) = A(n,ni) + a2(z);
    elseif mu == 3 %plasmid replication
           A(n,ni) = A(n,ni) + a3(z);
    elseif mu == 4 %plasmid dilution
           A(n,ni) = A(n,ni) + a4(z);
    end
end

%PROPENSITIES 

function out = a1(z)
global sysP
a = sysP.rcp(1);
out = a*z(2);
end

function out = a2(z)
global sysP
b = sysP.rcp(2);
out = b*z(1);
end

function out = a3(z)
global sysP
aP = sysP.rcp(4);
out = aP*z(2);
end

function out = a4(z)
global sysP
bP = sysP.rcp(5);
out = bP*z(2);
end

