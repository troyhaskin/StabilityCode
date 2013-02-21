function [] = UpdateThermodynamics(q)
    
    % Get pertinent data from the Database
    [Irho,Irhoi,Tguess] = Database('get','Irho','Irhoi','T');
    
    rho  = q(Irho)  ;
    rhoi = q(Irhoi) ;

    i = rhoi ./ rho;
    [T,PhaseCheck] = Temperature(rho,i,Tguess);
    
    Database('set','T'        , T                             ,...
                   'PhaseCheck', PhaseCheck                    ,...
                   'P'         , Pressure(rho,T,PhaseCheck)    ,...
                   'mu'        , Viscosity(rho,T,PhaseCheck)   );
end