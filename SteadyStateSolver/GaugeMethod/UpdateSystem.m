function [] = UpdateSystem(q)
    
    % Get pertinent data from the Database
    [Irho,Irhoi,Tguess,Ncv] = Database('get','Irho','Irhoi','T','Ncv');
    
    rho  = q(Irho)  ;
    rhoi = q(Irhoi) ;
    
    i = rhoi ./ rho;
    [T,PhaseCheck] = Temperature(rho,i,Tguess);
    P =  Pressure(rho,T,PhaseCheck);
    Database('set','T'        , T                           ,...
        'PhaseCheck', PhaseCheck                 ,...
        'P'         ,  P                         ,...
        'mu'        , Viscosity(rho,T,PhaseCheck));
    
    
    rhouNow = q(Ncv+1);
    [rhouOld,DeltaPold,DeltaPfric,dx] = Database('get','rhouNew','DeltaPold','DeltaPfric','dx');
    DeltaP = (P(end)-dx(end)/2*sum(DeltaPfric)) - P(1);
    
    switch(isnan(DeltaPold))
        case false
            
            rhouNew = rhouNow - DeltaP./(DeltaP-DeltaPold)*(rhouNow-rhouOld);
            
        case true
            
            % Half of double the speed for the first iteration
            if DeltaP > 0
                rhouNew = 1.1*rhouNow;
            else
                rhouNew   = 0.9*rhouNow;
            end
    end
    
    if rhouNew > rhouNow
        Outcome = 'increased';
    else
        Outcome = 'decreased';
    end
    
    % Store these values for later iterations
    rhouOld   = rhouNow;
    DeltaPold = DeltaP;
    
    Database('set','rhouNew',rhouNew,'rhouOld',rhouOld,'DeltaPold',DeltaPold);
    
    
    fprintf('Momentum was %s. DeltaP is now %+17.8E.\n',Outcome,DeltaP);
end