function DiffStateRate = DifferentiationFunction(x,parasDiffFunction)
    
    rH = parasDiffFunction(1);
    Kh = parasDiffFunction(2);
    nH = parasDiffFunction(3);
    
    %Hill function differentiation
    DiffStateRate = rH*(x.^nH./(x.^nH+Kh^nH));
    %DiffStateRate = x.^nH./(x.^nH+Kh^nH);
%     %linear differentiation
%     DiffStateRate = rH*(x+max(0,(x-55)*15))/Kh;
%     DiffStateRate = rH*x/Kh;
    
end











