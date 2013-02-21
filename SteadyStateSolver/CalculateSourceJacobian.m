function SourceJacobian = CalculateSourceJacobian(q,State,Info)
    
    % Allocate
    SourceJacobian = zeros(length(q));
    
    % Conserved variables
%     rho  = q(1);
    
    % Geometry info
    Grav = Info.g;

    % Total derivative w.r.t to the conserved variables
    [dPdxFric_rho,dPdxFric_rhov,dPdxFric_rhoi] = tauWallDerivatives(q,State,Info);
    
    % Source Jacobian formation
    SourceJacobian(1,:) = [          0            ,        0       ,      0         ];
    SourceJacobian(2,:) = [Grav-dPdxFric_rho  , -dPdxFric_rhov , -dPdxFric_rhoi     ];
    SourceJacobian(3,:) = [          0            ,        0       ,      0         ];

end