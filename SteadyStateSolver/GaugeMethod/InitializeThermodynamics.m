function [] = InitializeThermodynamics(q)
    
    % Get pertinent data from the Database
    [Irho,Irhoi] = Database('get','Irho','Irhoi');
    
    rho  = q(Irho)  ;
    rhoi = q(Irhoi) ;

    i = rhoi ./ rho;
    [T,PhaseCheck] = Temperature(rho,i);
    
    Database('insert','T'        , T                             ,...
                     'PhaseCheck', PhaseCheck                    ,...
                     'P'         , Pressure(rho,T,PhaseCheck)    ,...
                     'mu'        , Viscosity(rho,T,PhaseCheck)   );
end