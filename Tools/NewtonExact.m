function [x,varargout] = NewtonExact(x0,F,J)
    
    Tolerance = 1E-10   ;
    NotDone   = true    ;
    xk        = x0      ;
    
    while NotDone

        dx = J(xk) \ F(xk)  ;
        xk = xk - dx        ;
        
        NotDone = norm(dx,1) > Tolerance;

    end
    
    x = xk;
    
    if (nargout >= 2)
        varargout{1} = norm(F(xk),2);
    end
    
end