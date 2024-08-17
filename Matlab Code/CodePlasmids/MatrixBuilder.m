function [A] = MatrixBuilder(truncX,truncP,paras,uIn); 
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
rH = paras(7);
Kh = paras(8);
nH = paras(9);


s=(xmax+1)*(ymax+1)*2+1;
A=zeros(s);

%system parameters
sysP.rcp = [a  b burst aP bP pKill rH  Kh  nH  uIn];

sysP.S = [1 -1 0 0 0;  %s1
          0 0 1 -1 0
          0 0 0  0 1]; %s2
sysP.name = 'BS';
sysP.maxS = abs(max(max(sysP.S)));  %maximum S deviation

S = sysP.S;
M = size(sysP.S,2);


for n=1:s-1
        
     %compute coordinates of state n to align state space
     z = n2state2dPlus(n,xmax,s);
     
     %for off-diagonals
     for mu=1:M
         zj = z + S(:,mu);
         if zj(1)>=0 && zj(1)<=xmax && zj(2)>=0 && zj(2)<=ymax %guarding truncation boundary
             nj = state2n2dPlus(zj,xmax,s);
             A = addprop(z,n,nj,A,mu,xmax,s);
         end
     end

     %compute the diagonal (include J' values)   
     A(n,n) = -sum(A(n,:)); 
end

for n = 1:s-1
    %compute coordinates of state n
    z = n2state2dPlus(n,xmax,s);
    if z(2) == 0 %no plasmids
        OutRate = pKill;
    else
        OutRate = 0;
    end
    A(n,n) = A(n,n) - OutRate;
end

%crop the J' values
A(end,:)=[]; A(:,end)=[];

%transpose into general form
A = A.';

end