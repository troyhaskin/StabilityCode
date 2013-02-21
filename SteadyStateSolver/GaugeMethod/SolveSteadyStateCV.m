function [qNext,State] = SolveSteadyStateCV(qLast,Info,qGuess)
    
    if (nargin >= 3) && not(any(qGuess==0))
        qN = qGuess;
    else
        qN = 1.001*qLast;
    end
    
    [Flast,State]  = CalculateFlux(qLast,Info.Tguess);
    [Slast,State]  = CalculateSource(qLast,State,Info);
    
    NotDone = true;
    while NotDone
        [Residual,~,State] = CalculateResidual(qN,qLast,Flast,Slast,State,Info);
        ResJac    = CalculateResidualJacobian(qN,0,State,Info);
        ResJacInv = Invert3by3(ResJac);
        
        dqN = ResJacInv * Residual;
        qN  = qN - dqN;
        
        Error   = norm(dqN./qN,1) ;
        NotDone = Error > 1E-8;
        
    end
    
    State = UpdateState(qN,State,Info);
    qNext = qN;
    
end