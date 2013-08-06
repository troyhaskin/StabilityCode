function [x,Residual,LastResidualVec] = GMRESHouseholder(A,b,x0,Nrestarts,r0)
    
    % ============================================================= %
    %                             Set-Up                            % 
    % ============================================================= %
    
    % Length of all vectors
    N = length(b);
    
    %   Number of iterations before restarting calculation; i.e.,
    %   maximum number of update basis vectors to store before "forgeting" them
    %   and starting the GMRES process again with Nrestart's x0 and r0.
    if (nargin < 4) || isempty(Nrestarts)
        Nrestarts = N;
    end
    
    % Initial solution guess
    if (nargin < 3) || isempty(x0)
        x0 = zeros(N,1);
    end
    
    % Initial residual
    if (nargin < 5) || isempty(r0)
        r0 = b - A*x0;
    end
    
    SolveTolerance = 1E-10;
    
    % Matrix allocation
    Z = zeros(N,Nmax)   ; % Hold's update basis vectors
    H = zeros(N,Nmax)   ; % Holds Householder vectors for projections
    V = zeros(N,Nmax)   ; % Holds unitary matrix columns of QR decomposition
    R = zeros(N,Nmax)   ; % Holds upper-triangular matrix for least-square problem

    
    rNorm  = norm(r0,2)         ;
    Z(:,1) = r0 / rNorm         ;
    alpha  = r0*0               ;
    e      = [1 ; zeros(N-1,1)] ;
    
    % QR Factorization
    
    % Compute A*z1
    R(:,1) = A*Z(:,1)       ;
    
    % Compute Householder vector to bring R(:,1) into upper triangular form
    x      = R(:,1);
    u      = -sign(x(1)) * norm(x,2) * e(1:N) - x;
    H(:,1) = u / norm(u,2);
    
    % Apply projection to R to bring it into upper triangular form
    R(:,1) = R(:,1) - 2 * H(:,1) * (H(:,1)'*R(:,1));
    
    % Get the first column of the unitary matrix
    V(:,1) = e - 2 * H(:,1) * (H(:,1)'*e);
    
    % Residual update
    alpha(1) = V(:,1)'*r0           ;
    rk       = r0 - alpha(1)*V(:,1) ;
    
    rkm1Norm = rNorm        ;
    rkNorm   = norm(rk,2)   ;
    
    nu = 0.5;
    
    for k = 2:Nmax
        
        if rkNorm <= nu*rkm1Norm
            Z(:,k) = rk/rkNorm;
        else
            Z(:,k) = V(:,k-1);
        end
        
        % Compute and store A*zk in R
        R(:,k) = A*Z(:,k)   ;
        
        % Apply all previous projections to new the column
        for m = 1:k-1
            h      = H(1:N-m+1,m);
            R(:,k) = [R(1:m-1,k) ; R(m:N,k) - 2*h*(h'*R(m:N,k))];
        end
        
        % Get the next Householder vector
        x            = R(k:N,k)                                 ;
        h            = -sign(x(1)) * norm(x,2) * e(1:N-k+1) - x ;
        h            = h / norm(h)                              ;
        H(1:N-k+1,k) = h                                        ;
        
        %   Apply projection to R to bring it into upper triangular form;
        %   The triu() call explicitly zeros all strictly lower triangular
        %   components to minimize FP error.
        R(:,1:k) = triu([R(1:k-1,1:k) ; R(k:N,1:k) - 2 * h * (h'*R(k:N,1:k))]);
        
        % Get the k-th column of the current unitary matrix
        V(:,k) = [zeros(k-1,1) ; e(1:N-k+1) - 2*h*(h'*e(1:N-k+1))];
        for m = k-1:-1:1
            hm     = H(1:N-m+1,m);
            V(:,k) = [V(1:m-1,k) ; V(m:N,k) - 2*hm*(hm'*V(m:N,k))];
        end
        
        % Update residual
        alpha(k) = V(:,k)'*rk               ;
        rk       = rk - alpha(k)*V(:,k)     ;
        
        rkm1Norm = rkNorm;
        rkNorm   = norm(rk,2);
        
        % Solve least-squares problem
        if rkNorm < SolveTolerance
            break;
        end
        
    end
    
    % Solve the least-squares problem
    yk = R(1:k,1:k) \ alpha(1:k);
    
    % Calculate actual shift of guess
    dx = Z(:,1:k) * yk  ;
    x  = x0 + dx;

end
