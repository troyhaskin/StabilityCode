function FluxJacobian = CalculateFluxJacobian(q,State)
    
    % Allocate
    FluxJacobian = zeros(length(q));
    
    % Conserved variables
    rho  = q(1);
    rhov = q(2);
    rhoi = q(3);
    
    % Primitives
    v = rhov / rho;
    Int = rhoi / rho;
    
    % Thermodynamics
    T = State.T;
    P = State.P;
    h = Int + P/rho;

    % Total derivative w.r.t to the conserved variables
    [dPdrho,dPdi] = Pressure_TotalDerivatives(rho,T,Int,State.PhaseCheck);
    
    % Flux Jacobian formation
    FluxJacobian(1,:) = [        0      ,  1  ,      0    ];
    FluxJacobian(2,:) = [ dPdrho - v^2  , 2*v ,    dPdi   ];
    FluxJacobian(3,:) = [v*(dPdrho - h) ,  h  , v*(1+dPdi)];

end