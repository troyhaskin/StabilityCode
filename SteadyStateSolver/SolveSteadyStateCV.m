function [qNext,State] = SolveSteadyStateCV(qLast,Info,qGuess)
    
    if (nargin >= 3) && not(any(qGuess==0))
        qN = qGuess;
    else
        qN = 1.001*qLast;
    end
    
    % Pre-calculate fixed Flux and Source values from previous iteration
    [Flast,State]  = CalculateFlux(qLast,Info.Tguess,Info.PhaseCheck);
    [Slast,State]  = CalculateSource(qLast,State,Info);
    
    % Calculate initial search direction and error
    [Residual,~,State] = CalculateResidual(qN,qLast,Flast,Slast,State,Info) ;
    ResJac             = CalculateResidualJacobian(qN,0,State,Info)         ;
    ResJacInv          = Invert3by3(ResJac)                                 ;
    dqN                = ResJacInv * Residual                               ;
    NormResidual       = norm(Residual,1)                                   ;
    
    % Loop conditional
    NotDone = true;
    
    while NotDone
        
        qCandidate     = qN - dqN  ;
        alpha          = 1         ;
        StepSizeTooBig = true      ;
        
        while StepSizeTooBig
            [Residual,~,State]    = CalculateResidual(qCandidate,qLast,Flast,Slast,State,Info);
            NormResidualCandidate = norm(Residual,1);

            ErrorNotReduced = (NormResidualCandidate > NormResidual);
            StepSizeTooBig  = ErrorNotReduced || isnan(NormResidualCandidate);

            if StepSizeTooBig
                alpha      = alpha / 2;
                qCandidate = qN - alpha * dqN  ;
                disp(['alpha reduced to :',num2str(alpha,'%6E')]);
            end

            if alpha < eps();
                break;
            end;

        end
        qN           = qCandidate               ;
        NormStepSize = norm(alpha*dqN./qN,1)    ;
        NormResidual = NormResidualCandidate    ;

        ResJac    = CalculateResidualJacobian(qN,0,State,Info);
        ResJacInv = Invert3by3(ResJac);
        dqN = ResJacInv * Residual; % Search direction
        
        NotDone = NormStepSize > 1E-8;
        
    end
    
    State = UpdateState(qN,State,Info);
    qNext = qN;
    
end