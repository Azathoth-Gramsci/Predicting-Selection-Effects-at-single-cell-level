function [A] = MatrixBuilderNoDiff(truncX,truncP,paras); 
global s
global sysP
global xmax

xmax = truncX;
ymax = truncP;

a = paras(1);
b = paras(2);
burst = paras(3);
aP = paras(4);
bP = paras(5);
pKill = paras(6);

s=(xmax+1)*(ymax+1)+1;
A=zeros(s);

%system parameters
sysP.rcp = [a  b burst aP bP pKill];

sysP.S = [1 -1 0 0;  %s1
          0 0 1 -1]; %s2
sysP.name = 'BS';
sysP.maxS = abs(max(max(sysP.S)));  %maximum S deviation

S = sysP.S;
M = size(sysP.S,2);

for n=1:s-1 
     %compute coordinates of state n to align state space
     z = n2state(n,xmax);     
     for mu=1:M
         zj = z + S(:,mu);
         if zj(1)>=0 && zj(1)<=xmax && zj(2)>=0 && zj(2)<=ymax %guarding truncation boundary
             nj = state2n(zj,xmax,ymax);
             A = addpropNoDiff(z,n,nj,A,mu,xmax,ymax);
         end
     end
     %compute the diagonal
     A(n,n) = -sum(A(n,:)); 
end

for n = 1:s-1
    %compute coordinates of state n
    z = n2state(n,xmax);
    if z(2) == 0 %no plasmids
        OutRate = pKill;
    else
        OutRate = 0;
    end    
    A(n,n) = A(n,n) - OutRate;
end

%crop values
A(end,:)=[]; A(:,end)=[];

%transpose into general form
A = A.';

end