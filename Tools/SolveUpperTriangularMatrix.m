function x = SolveUpperTriangularMatrix(U,b)

    N = length(b);
    x = b;

    x(N) = b(N) / U(N,N);
    
    for k = (N-1):-1:1
        x(k) = (b(k) - U(k,k+1:N)*x(k+1:N))/U(k,k);
    end
    
    
end