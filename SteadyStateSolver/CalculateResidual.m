function [Residual,qMid,State]  = CalculateResidual(qN,qLast,Flast,Slast,State,Info)
    
    
    % Form the sources at the CV's center
%     qMid       = (qN+qLast)/2;
%     [Smid,~] = CalculateSource(qMid,State,Info);
    qMid = 0;
    
    % Form the fluxes at CV's boundaries
    [Fn,State] = CalculateFlux(qN,State.T,State.PhaseCheck);
    % Form the sources at the CV's boundaries
    [Sn,State] = CalculateSource(qN,State,Info);
    
    % Form Residual
    Residual = ((Fn - Flast) - Info.dx/2*(Sn + Slast));
%      Residual = (Fn - Flast) - Info.dx/6*(Sn + 4*Smid + Slast);
    
end