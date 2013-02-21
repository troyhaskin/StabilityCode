function R = CalculateResidual(q)
    
    
%     % Form the sources at the CV's center
%     qMid = 0;
    
    [MaskL,MaskR,dx] = Database('get','MaskL','MaskR','dx');
    
    % Flux values
    FluxAll = HEMFlux(q)        ;
    FluxL   = FluxAll(MaskL)    ;
    FluxR   = FluxAll(MaskR)    ;
    
    % Source values
    SourceAll = HEMSource(q)        ;
    SourceL   = SourceAll(MaskL)    ;
    SourceR   = SourceAll(MaskR)    ;
    
    % Form Residual
    R = (FluxR - FluxL) - dx.*(SourceL + SourceR)/2;
    
end