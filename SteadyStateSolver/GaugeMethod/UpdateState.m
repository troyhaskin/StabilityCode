function State = UpdateState(q,State,Info)
   
    % Conserved variables
    rho  = q(1);
    rhoi = q(3);
    
    % Primitives
    Int = rhoi/rho;
    
    % Thermodynamics
    State.T = Temperature(rho,Int,State.T,false);
    State.P = Pressure   (rho,State.T,false);
    
    % Transport
    [~,State] = tauWall(q,State,Info);
    
end