function x = LowerBidiagonalSolve(A,b)
    x = b;
    x(1) = b(1)/A(1,1);
    
    N = length(b);
    
    for k = 2:N
        x(k) = (b(k) - A(k,k-1)*x(k-1))/A(k,k);
    end

end