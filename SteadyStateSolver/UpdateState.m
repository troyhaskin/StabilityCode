function State = UpdateState(q,State,Info)
   
    % Conserved variables
    rho  = q(1);
    rhoi = q(3);
    
    % Primitives
    Int = rhoi/rho;
    
    % Thermodynamics
    State.T = Temperature(rho,Int,State.T);
    State.P = Pressure   (rho,State.T);
    
    % Transport
    [~,State] = tauWall(q,State,Info);
    
end