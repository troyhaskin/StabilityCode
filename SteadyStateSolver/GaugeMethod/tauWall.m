function [DeltaP,State] = tauWall(q,State,Info)
    
    % Conserved variables
    rho  = q(1);
    rhov = q(2);
    
    % Primitives
    v = rhov./rho;
    
    % Geometry
    Lchar  = Info.Lchar;
    EpsRel = Info.EpsRel;
    
    % Thermodynamics properties
    T  = State.T;
    mu = Viscosity(rho,T,false);
    
    % Friction values
    Re       = ReynoldsNumber(rhov,mu,Lchar)          ;
%     [f,f_Re] = FrictionFactorCW(Re,EpsRel,Info.fGuess)            ;
    [f,f_Re] = FrictionFactorConstant(0.0137)          ;
    
    % Update the State information
    State.f    = f      ;
    State.Re   = Re     ;
    State.f_Re = f_Re   ;
    State.mu   = mu     ;
    
    % DeltaP due to friction
    DeltaP = 2 * f .* (1./Lchar) .* rho.*v.^2   ;
    
end