function [x,IterationsNonlinear] = JFNKHouseholder(x0,F,epsilon)

    
    % ================================================================= %
    %                               Set-Up                              %
    % ================================================================= %
    
    % Length
    N    = length(x0)   ;
    Nmax = N            ;
    
    % Tolerances
    LinearTolerance    = 1E-10;
    NonlinearTolerance = 1E-10;

    % Matrix allocation
    Z = zeros(N,Nmax)   ; % Hold's update basis vectors
    H = zeros(N,Nmax)   ; % Holds Householder vectors for projections
    V = zeros(N,Nmax)   ; % Holds unitary matrix columns of QR decomposition
    R = zeros(N,Nmax)   ; % Holds upper-triangular matrix for least-squares problem
    
    % Unit vector used for projections
    e     = [1 ; zeros(N-1,1)]  ;
    alpha = zeros(N,1)          ;

    % Threshold parameter used to determine which basis to use in the GMRES iterations
    nu = 0.85;




    % ================================================================= %
    %                            JFNK Iteration                         %
    % ================================================================= %

    % Initial r0
    r0    = -F(x0)      ;
    rNorm = norm(r0,2)  ;
    
    % Backtrack relaxor
    relaxor = 0.5;
    
    % Determine if the loop is needed
    NotDone = rNorm > NonlinearTolerance;
    
    % Let x = x0
    x  = x0;
    r  = r0;
    
    % Counters
    IterationsNonlinear = 0;
    
    while NotDone
        
        % Solve linear system
        Nstop = GMRES(x,r,rNorm);   % Solve the linear system to LinearTolerance
        I     = 1:Nstop         ;   % Vector of Indices for the solve
        
        % Update x
        yk = R(I,I) \ alpha(I)      ;   % Solve the least-squares problem
        dx = Z(:,I) * yk            ;   % Calculate full Newton update
        
        % Backtracker
        while (norm(F(x),2) < norm(F(x + dx),2))
            dx = relaxor * dx;
        end
        x  = x + dx ;   % Calculate relaxed x value
        
        % Check non-linear residual
        r     = -F(x)       ;
        rNorm = norm(r,2)   ;

        Show( rNorm );

        % Loop break check
        NotDone = rNorm > NonlinearTolerance;

        IterationsNonlinear = IterationsNonlinear + 1;
    end
    
    
    
    
    % ================================================================= %
    %                          GMRES SubFunction                        %
    % ================================================================= %
    function Nstop = GMRES(xk,rk0,rk0Norm)
        
        Z(:,1) = rk0 / rk0Norm ; % First basis vector for update
        
        % First Step (k = 1)
        % Compute J*z1 and store in R
        R(:,1) = (F(xk + epsilon*Z(:,1)) + rk0) / epsilon;
        
        % Compute Householder vector to bring R(:,1) into upper triangular form
        h      = R(:,1);
        h      = -sign(h(1)) * norm(h,2) * e(1:N) - h;
        H(:,1) = h / norm(h,2);
        
        % Apply projection to R to bring it into upper triangular form
        R(:,1) = R(:,1) - 2 * H(:,1) * (H(:,1)'*R(:,1));
        
        % Get the first column of the unitary matrix
        V(:,1) = e - 2 * H(:,1) * (H(:,1)'*e);
        
        % Residual update
        alpha(1) = V(:,1)'*rk0           ;
        rk       = rk0 - alpha(1)*V(:,1) ;
        
        % Assign residual norms to determine which basis to use
        rkm1Norm = rk0Norm      ;
        rkNorm   = norm(rk,2)   ;
        
        
        for k = 2:Nmax
            
            % Choose the next basis vector
            if rkNorm <= nu*rkm1Norm
                Z(:,k) = rk/rkNorm;     % GCR (RB-SGMRES) basis
            else
                Z(:,k) = V(:,k-1);      % Simpler GMRES basis
            end
            
            % Compute and store A*zk in R
            R(:,k) = (F(xk + epsilon*Z(:,k)) + rk0) / epsilon;
            
            % Apply all previous projections to new the column
            for m = 1:k-1
                h      = H(1:N-m+1,m);
                R(:,k) = [R(1:m-1,k) ; R(m:N,k) - 2*h*(h'*R(m:N,k))];
            end
            
            % Get the next Householder vector
            h            = R(k:N,k)                                 ;
            h            = -sign(h(1)) * norm(h,2) * e(1:N-k+1) - h ;
            h            = h / (norm(h) + eps(h))                   ;
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
            
            % Update residual norms
            rkm1Norm = rkNorm;
            rkNorm   = norm(rk,2);
            
            % Solve least-squares problem
            if rkNorm < LinearTolerance
                break;
            end
            
        end
        
        Nstop    = k        ;
    end
end
