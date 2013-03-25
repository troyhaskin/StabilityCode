function [F,State] = CalculateFlux(q,Tguess,PhaseCheck)
    
    if (nargin < 2)
        Tguess = [];
    end
    
    
    % Assign conserved variables
    rho  = q(1);
    rhov = q(2);
    rhoi = q(3);
    
    % Construct primitives
    v = rhov / rho;
    i = rhoi / rho;
    
    % Theromdynamics properties
    [T,PhaseCheck] = Temperature(rho,i,Tguess);
%     Show([rho,i,T]','%+17.8E');
    P = Pressure(rho,T,PhaseCheck);
    
    % Form flux vector
    F = [      rhov         ;...
        v *  rhov + P     ;...
        v * (rhoi + P)    ];
    
    % Form state struct
    State.P = P;
    State.T = T;
    State.PhaseCheck = PhaseCheck;
    
end