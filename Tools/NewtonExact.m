function [x,varargout] = NewtonExact(x0,F,J)
    
    Tolerance  = 1E-10  ;
    NotDone    = true   ;
    xk         = x0     ;
    Iterations = 0      ;
    
    while NotDone

        dx = J(xk) \ F(xk)  ;
        xk = xk - dx        ;
        
        Show(dx);
        
        NotDone    = norm(F(xk),2) > Tolerance ;
        Iterations = Iterations + 1         ;
    end
    
    x = xk;
    
    if (nargout >= 2)
        varargout{1} = norm(F(xk),2);
    end
    
    if (nargout >= 3)
        varargout{2} = Iterations;
    end
    
end