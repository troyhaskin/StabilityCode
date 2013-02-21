function S_q = CalculateSourceJacobian(q)
    
    % Database call
    Requests = {'Irho','Irhou','Irhoi','T','P','PhaseCheck'};
    [Irho,Irhou] = Database('get',Requests{:});
    
    % Database call
    Requests = 'S_q';
    S_q = Database('get',Requests{:});
    
    % Conserved variables
    rho  = q(Irho);

    % Total derivative w.r.t to the conserved variables
    [dPdxFric_rho,dPdxFric_rhov,dPdxFric_rhoi] = tauWallDerivatives(q);
    
    % Source Jacobian formation
    S_q(Irhou,1:3:end) = Gravity.*rho-dPdxFric_rho ; -dPdxFric_rhov ; -dPdxFric_rhoi ];

end