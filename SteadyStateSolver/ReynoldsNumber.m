function Re = ReynoldsNumber(rhov,mu,Lchar)
    
    Re = rhov .* Lchar ./ mu;
    
end