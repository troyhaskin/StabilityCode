function [x,Residual,LastResidualVec] = GMRESloops(A,b,x0,Nmax,r0)
    
    N = length(b);
    
    if (nargin < 4) || isempty(Nmax)
        Nmax = length(b);
    end
    
    if (nargin < 3) || isempty(x0)
        x0 = zeros(length(b),1);
    end

    if (nargin < 5) || isempty(r0)
        r0 = b - A*x0;
    end

    Nallocate = (Nmax == N)*N + (Nmax < N)*(Nmax+1);

    beta     = norm(r0,2) * [1;zeros(Nallocate,1)]  ;

    H        = zeros(Nallocate , Nmax)  ;
    Q        = zeros(N         , Nmax)  ;
    Q(:,1)   = r0./norm(r0,2)           ;

    Residual = zeros(Nmax,1);
    
    Niter = Nmax;

    for k = 1:Niter

        % Form 
        v = A*Q(:,k);
        
        % Calculate the column of the Hessenberg for the current Unitary set of vecotrs
        for m = 1:k
            H(m,k) = v'*Q(:,m);
        end
        
        % Project H*Q out of the current v
        for m = 1:k
            v = v - H(m,k)*Q(:,m);
        end

        H(k+1,k) = norm(v,2)     ;
        Q(:,k+1) = v / H(k+1,k)  ;

        % Create QR decomposition of current Hessenberg matrix
        [Omega,R] = qr(H(1:k+1,1:k),0);
        
        % Solve least-squares problem
        yk = R \ ( Omega' * beta(1:k+1) );
        
        % Compute iterate's update
        dx = Q(:,1:k)*yk;
        
        % Compute the new residual
        rk     = r0 - A*dx;
        
        % Compute and store the residual
        rkNorm      = norm(rk,2)    ;
        Residual(k) = rkNorm        ;
    end
    
    Residual = [norm(r0,2);Residual];
    x =  (x0 + dx);
    LastResidualVec = rk;
% %     if PerformFinalProjection
%         H(:,Nmax)   = Q'*v                 ;
%         [Omega,R]   = qr(H,0)    ;
%         yk          = R\(Omega'*beta)           ;
%         dx          = Q*yk              ;
%         rk          = r0 - A*dx         ;
%         rkNorm      = norm(rk,2)        ;
%         Residual(Nmax) = rkNorm/norm(r0,2) ;
% %     end
%     
%     x = x0 + dx;
    
end
