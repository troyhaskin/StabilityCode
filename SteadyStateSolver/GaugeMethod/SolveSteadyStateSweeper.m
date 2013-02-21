function [q,Thermo] = SolveSteadyStateSweeper(q,Thermo,Info)
    
%    Buoyancy = 0;
%    Friction = 0;
    
    % Sweep through the control volumes
    for k = 1 : Info.Ncv
    
        % Get the next control volume's state
        [q(:,k+1),State] = SolveSteadyStateCV(q(:,k),Info.CV(k),q(:,k+1));
        
        % Push properties
        Thermo.rho(k+1) = q(1,k+1)          ;
        Thermo.i(k+1)   = q(3,k+1)/q(1,k+1) ;
        Thermo.P(k+1)   = State.P           ;
        Thermo.T(k+1)   = State.T           ;
        Thermo.mu(k+1)  = State.mu          ;
        Thermo.Re(k+1)  = State.Re          ;
        Thermo.f(k+1)   = State.f           ;
        
        % Update Temperature guess value
        Info.CV(k+1).Tguess = State.T;
        Info.CV(k+1).fGuess = State.f;
        
%        % Perform Buoyancy and  Friction force calculation
%        Lchar = Info.CV(k).Lchar;
%        dx    = Info.CV(k).dx;
%        
%        Frictionk = 2*Thermo.f(k)*q(2,k)^2.*(1/Lchar/Thermo.rho(k));
%        Buoyancyk = Info.CV(k).g * Thermo.rho(k);
%        
%        Frictionkp1 = 2*Thermo.f(k+1)*q(2,k+1)^2.*(1/Lchar/Thermo.rho(k+1));
%        Buoyancykp1 = Info.CV(k).g * Thermo.rho(k+1);
%        
%        Buoyancy = Buoyancy + dx/2 * (Frictionk + Frictionkp1);
%        Friction = Friction + dx/2 * (Buoyancyk + Buoyancykp1);
    end
    
%    % Net Buoyancy and Friction forces
%    Thermo.Buoyancy = Buoyancy;
%    Thermo.Friction = Friction;
end




