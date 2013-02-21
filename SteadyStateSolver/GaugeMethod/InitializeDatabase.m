function [] = InitializeDatabase(q)
    
    [Sizer,Ncv] = Database('get','SizerLaw','Ncv');
    
    Database('insert','f',1E-4.*Sizer,'Re',Sizer,'dfdRe',Sizer,'rhouOld',NaN,'rhouNew',q(Ncv+1),'DeltaPold',NaN,'DeltaPfric',0);
    InitializeThermodynamics(q);

end


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