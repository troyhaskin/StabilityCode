function ResidualJacobian = CalculateResidualJacobian(q,qMid,State,Info)
    
    % CV spacing
    dx = Info.dx;
    
    % Form all the Jacobians
    FluxJacobian      = CalculateFluxJacobian  (q   ,State);
    SourceJacobian    = CalculateSourceJacobian(q   ,State,Info);
%     SourceJacobianMid = CalculateSourceJacobian(qMid,State,Info);
    
    % Form the residual Jacobian
%     ResidualJacobian = FluxJacobian - dx/6*(2*SourceJacobianMid + SourceJacobian);
    ResidualJacobian = FluxJacobian - dx/2 * SourceJacobian;
    
end