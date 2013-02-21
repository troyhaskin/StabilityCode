function [fSol,varargout] = FrictionFactorCW(Reynolds,EpsRelative,fGuess)
    
    Rez   = @(fVar,Re,EpsRel) CWResidual(EpsRel,Re,fVar);
    DRez  = @(fVar,Re,EpsRel) CWResidualFirstDerivative(EpsRel,Re,fVar);
    
    if (nargin < 3) || isempty(fGuess)
        f = 1E-3 + 0 * Reynolds      ;
    else
        if not(isscalar(fGuess))
            f = fGuess;
        else
            f = fGuess + 0 * Reynolds ;
        end
    end
    
    fSol            = f                        ;
    ErrL2Rel        = 1    + 0 * Reynolds      ;
    ErrNewtonStop   = 100 * eps()              ;
    Iter            = 0                        ;
    IterStop        = 1E4                      ;
    Mask            = ErrL2Rel > ErrNewtonStop ;
    NotConverged    = true                     ;
    
    
% ====================================================================== %
%                          Newton Iterations                             %
% ====================================================================== %
    while NotConverged
        RezIter  = @(fVar) Rez(fVar,Reynolds(Mask),EpsRelative(Mask));
        DRezIter = @(fVar) DRez(fVar,Reynolds(Mask),EpsRelative(Mask));
        
        fSol(Mask) = NewtonUpdate(f(Mask),RezIter,DRezIter);
        
        ErrL2Rel = sqrt((f-fSol).^2)./abs(f);
        Iter     = Iter + 1                 ;
        Mask     = ErrL2Rel > ErrNewtonStop ;
        
        TooManyIter      = Iter > IterStop                  ;
        NotAllfConverged = any(Mask)                        ;
        NotConverged     = NotAllfConverged || TooManyIter  ;
        
        f	= fSol;
    end
% ====================================================================== %
    



% ====================================================================== %
%                        Friction Factor Derivative                      %
% ====================================================================== %
% If two output arguments are expected, the derivative of 
% the friction factor with respect to Reynolds Number is returned.
%
% This derivative can be obtained explicitly (using the known friction 
% factors) from the Colebrook-White equation via Chain Rule.

    if(nargout == 2)
        [~,C2,C3,C4] = GetCoefficients();
        
        fSqrt = sqrt(fSol);
        Numer = 2 * C2 * C4 * fSol.^(3/2);
        Denom = Reynolds.*(C2*C4*fSqrt + log(10)*(C4 + C3*EpsRelative.*Reynolds.*fSqrt));
        
        varargout{1} = -Numer./Denom;
    end
    
end


% ============================================================================ %
%                    Required Functions for Newton Iteration                   %
% ============================================================================ %
% ---------- Coefficients used in Colebrook-White Equation ---------- %
function [C1,C2,C3,C4] = GetCoefficients()
	C1  = 3.48;
    C2  = 4.00;
    C3  = 2.00;
    C4  = 9.35;
end


% ---------- Colebrook-White Equation in Residual Form ---------- %
function Residual = CWResidual(EoD,Re,f)

    [C1,C2,C3,C4] = GetCoefficients();
    
    fInvSqrt    = 1./(sqrt(f))                                              ;
    Residual    = fInvSqrt - (C1 - C2 * log10(C3*EoD + C4./Re.*fInvSqrt))	;
    
end

% ---------- Derivative of Colebrook-White Equation in Residual Form  with respect to f ---------- %
function Residual = CWResidualFirstDerivative(EoD,Re,f)
    
    [~,C2,C3,C4] = GetCoefficients();
    
    fTo3half    = f.^(3/2);
    Ln100       = log(100);
    Residual    = -(0.5 ./fTo3half + C2*C4./(C4*f*Ln100 + C3*EoD.*fTo3half.*Re*Ln100));
end

% ---------- Newton's Method ---------- %
function xNext = NewtonUpdate(xn,f,Df)
    
    xNext = xn - f(xn)./Df(xn);
    
end


% ============================================================================ %
%                    Required Functions for Halley Iteration (if wanted)       %
% ============================================================================ %
% ---------- Second derivative of Colebrook-White Equation in Residual Form  with respect to f ---------- %
function Residual = CWResidualSecondDerivative(EoD,Re,f) %#ok<DEFNU>
    
    [~,C2,C3,C4] = GetCoefficients();
    
    fSqrt       = sqrt(f);
    Ln10        = log(10);
    Term1       = C2*(2*C4^2.*fSqrt + 3*C3*C4*EoD*f.*Re);
    Term2       = (C4 + C3*EoD*fSqrt.*Re).^2*Ln10;
    Residual    = (0.75 + 0.25*Term1/Term2)./fSqrt.^5;
    
end


% ---------- Halley's Method ---------- %
function xNext = HalleyUpdate(xn,f,Df,D2f) %#ok<DEFNU>
    
    Update  = 2*f(xn).*Df(xn)./(2*Df(xn).^2-f(xn).*D2f(xn));
    xNext   = xn - Update;
    
end
