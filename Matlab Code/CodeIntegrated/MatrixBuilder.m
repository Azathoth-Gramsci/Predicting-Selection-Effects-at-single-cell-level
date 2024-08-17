function [A] = MatrixBuilder(truncX,paras,uIn); 
global s
global sysP
global xmax

xmax = truncX;

a = paras(1);
b = paras(2);
burst = paras(3);
rH = paras(4);
Kh = paras(5);
nH = paras(6);


s=(xmax+1)*2+1;
A=zeros(s);

%system parameters
sysP.rcp = [a  b burst rH  Kh  nH  uIn];

sysP.S = [1 -1 0;  %s1
          0 0 1]; %s2
sysP.name = 'BS';
sysP.maxS = abs(max(max(sysP.S)));  %maximum S deviation

S = sysP.S;
M = size(sysP.S,2);


for n=1:s-1
        
     %compute coordinates of state n
     z = n2state(n,xmax);
     %for off-diagonals
     for mu=1:M
         zj = z + S(:,mu);
         if zj(1)>=0 && zj(1)<=xmax && zj(2)>=0 && zj(2)<=1 %guarding truncation boundary
             nj = state2n(zj,xmax,2);
             A = addprop(z,n,nj,A,mu,xmax,s);
         end
     end

     %compute the diagonal (include J' values)   
     A(n,n) = -sum(A(n,:)); 
end

%crop the J' values
A(end,:)=[]; A(:,end)=[];

%transpose into general form
A = A.';

end