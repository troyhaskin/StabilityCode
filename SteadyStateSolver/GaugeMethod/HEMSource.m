function S = HEMSource(q)
    
    [Irho,Irhou,Irhoi,Gravity,qAdd] = Database('get','Irho','Irhou','Irhoi','Gravity','qAdd');
    
    rho  = q(Irho);
    rhou = q(Irhou);
    
    dPdxFric = tauWall(rho,rhou);
    
    S        =           q * 0              ;
    S(Irho)  =             0                ;
    S(Irhou) = Gravity .* rho - dPdxFric    ;
    S(Irhoi) =          qAdd                ;

end

function DeltaP = tauWall(rho,rhov)
    
    [mu,Lchar,EpsRel,fGuess] = Database('get','mu','Lchar','EpsRel','f');
    
    % Friction values
    Re = ReynoldsNumber(rhov,mu,Lchar)    ;
    [f,dfdRe] = FrictionFactorCW(Re,EpsRel,fGuess)      ;
    
    % DeltaP due to friction
    DeltaP = 2 * f .* (1./Lchar) .* rhov .* abs(rhov) ./ rho   ;
    
    Database('set','dfdRe',dfdRe,'f',f,'Re',Re,'DeltaPfric',DeltaP([1,end]));
end

