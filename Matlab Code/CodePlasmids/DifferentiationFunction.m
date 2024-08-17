function DiffStateRate = DifferentiationFunction(x,parasDiffFunction)
    
    rH = parasDiffFunction(1);
    Kh = parasDiffFunction(2);
    nH = parasDiffFunction(3);
    
    %Hill function differentiation
    DiffStateRate = rH*(x.^nH./(x.^nH+Kh^nH));
    
end











