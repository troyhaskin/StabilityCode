function Qnp1 = EvolveInTime(Q)
    
    % First-order method (time and space).
    dQ = GetTimeStepUpdates(Q);
    Qnp1 = Q - dQ;
    
end