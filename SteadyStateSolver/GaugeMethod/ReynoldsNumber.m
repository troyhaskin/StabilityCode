function Re = ReynoldsNumber(rhov,mu,Lchar)
    
    Re = abs(rhov) .* Lchar ./ mu;
    
end