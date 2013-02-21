function x = BiconjugateGradientMethod(A,x0,b)
    
    At = A';
    
    x  = x0;
    xs = x0;
    
    rk  = b - A*x0;
    rsk = rk;
    
    p  = rk;
    ps = rsk;
    
    for k = 1:length(x0)
        alpha = (rsk'*rk)/(ps'*A*p);
        
        x  = x  + alpha*p ;
        xs = xs + alpha*ps;
        
        rkp1  = rk  - alpha * A  * p ;
        rskp1 = rsk - alpha * At * ps;
        
        beta = (rskp1'*rkp1)/(rsk'*rk);
        
        p  = rkp1  + beta*p ;
        ps = rskp1 + beta*ps;
        
        rk  = rkp1 ;
        rsk = rskp1;
    end
    
end