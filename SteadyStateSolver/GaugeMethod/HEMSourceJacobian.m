function S_q = HEMSourceJacobian(q)
    
    % Database call
    Requests = {'S_q','LdF1dq','LdF2dq','LdF3dq','Gravity'};
    [S_q,LdS1dq,LdS2dq,LdS3dq,Gravity] = Database('get',Requests{:});
    
    % Conserved variables
    %     rho  = q(Irho);
    
    % Total derivative w.r.t to the conserved variables
    [dPdxFric_rho,dPdxFric_rhov,dPdxFric_rhoi] = tauWallDerivatives(q);
    
    % Source Jacobian formation
    S_q(LdS1dq) = 0;
    S_q(LdS2dq) = [Gravity - dPdxFric_rho ; -dPdxFric_rhov ; -dPdxFric_rhoi ];
    S_q(LdS3dq) = 0;
end