function Qnew = SolveImplicitSystem(Qold,dt)
    
    Qnew      = Qold  ;

    IterMax   = 1E2 ;
    Tolerance = 1E-6;

    Iter      = 0;
    NotDone   = true;
    
    while NotDone
        
        [Rez,RezJac] = SystemResidual(Qnew,Qold,dt); 
        
        dQ   = -RezJac \ Rez;
        Qnew = Qnew + dQ;
        
        Iter  = Iter + 1;
        Error = norm(dQ,2);
        
        NotConverged = Error > Tolerance            ;
        BelowIterMax = Iter  < IterMax              ;
        NotDone      = NotConverged & BelowIterMax  ;
        
    end
    
end