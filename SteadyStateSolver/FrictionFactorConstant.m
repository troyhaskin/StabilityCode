function [fSol,varargout] = FrictionFactorConstant(C,N)
    
    if (nargin < 2) || isempty(N)
        N = 1;
    end
    
    fSol = C * ones(N,1);
    
    if(nargout == 2)
        varargout{1} =fSol*0;
    end
    
end
