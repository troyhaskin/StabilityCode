function F_q = HEMFluxJacobian(q)
    
    % Database call
    Requests = {'Irho','Irhou','Irhoi','T','P','PhaseCheck'};
    [Irho,Irhou,Irhoi,T,P,PhaseCheck] = Database('get',Requests{:});
    
    % Database call
    Requests = {'F_q','LdF1drho','LdF1drhou','LdF1drhoi','LdF2dq','LdF3dq'};
    [F_q,LdF1drho,LdF1drhou,LdF1drhoi,LdF2dq,LdF3dq] = Database('get',Requests{:});

    % Conserved variables
    rho  = q(Irho);
    rhou = q(Irhou);
    rhoi = q(Irhoi);

    % Primitives
    u = rhou ./ rho;
    i = rhoi ./ rho;

    % Thermodynamics
    h = i + P./rho;

    % Total derivative w.r.t to the conserved variables
    [dPdrho,dPdi] = Pressure_TotalDerivatives(rho,T,i,PhaseCheck);

    % dF1dq
    F_q(LdF1drho)  =                     0                     ;
    F_q(LdF1drhou) =                     1                     ;
    F_q(LdF1drhoi) =                     0                     ;
    % dF2dq
    F_q(LdF2dq)    = [dPdrho - u.^2   ; 2*u ;     dPdi      ]  ;
    % dF3dq
    F_q(LdF3dq)    = [u.*(dPdrho - h) ;  h  ;  u.*(1+dPdi)  ]  ;

end


