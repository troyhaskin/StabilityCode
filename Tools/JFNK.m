function qk = JFNK(q0,F,epsilon)
    
    qk = q0         ;
    N  = length(qk) ;
    
    % Initial residual is how far
%     epsilon = 1E-6      ;
    r0      = F(qk)     ;
    Fqk     = r0        ;
    rNorm   = norm(r0,2);
    
    % GMRES arrays
    beta     = rNorm * [1;zeros(N,1)]   ; % Residual to be minimized
    H        = zeros(N+1, N )           ; % Iterative upper Hessenberg matrix allocation
    Q        = zeros( N ,N+1)           ; % Iterative Unitary          matrix allocation
    Q(:,1)   = r0 ./ rNorm              ; % Initial vector for the Unitary matrix

    % Iteration preparation
    NotDone   = true    ; % Loop breaker
    Tolerance = 1E-10   ; % Residual norm convergence criterion
    rBest     = rNorm   ; % Best residual norm acheived (used for backtracking)
    
    % Iteration counters
    IterationsNonlinear = 0;
    IterationsLinear    = 0;
    
    % Start: Nonlinear Solve
    while NotDone
        
        % Start: GMRES Linear Solve
        for k = 1:N
            
            % Start: Hessenberg Decomposition

                % Form the matrix-vector product A*Q_k via a finite difference
                v = (F(qk + epsilon*Q(:,k)) - Fqk) / epsilon;

                % Compute the Arnoldi vectors and associated Hessenberg columns by 
                % recurrence
                H(1:k,k) = Q(:,1:k)'*v; 
                v        = v - Q(:,1:k)*H(1:k,k);
                H(k+1,k) = norm(v,2);
                Q(:,k+1) = v/H(k+1,k);
            
            % End: Hessenberg Decomposition


            % Start: Least-Squares Minimization

                % Solve the least-squares minimization problem: min ||beta - H_k * y||
                [Omega,R] = qr(H(1:k+1,1:k),0);
                yk = R \ (Omega'*beta(1:k+1));

            % End: Least-Squares Minimization


            % Calculate the dq that will decrease the residual in the linear problem
            dx = Q(:,1:k)*yk;
            
            IterationsLinear = IterationsLinear + 1;
        end

        % Compute the new residual after making a full (quasi-) Newton step
        rNew = norm(F(qk - dx),2);

        % If the residual is worse, backtrack.
        if rNew > rBest

            alpha = 1;
            while rNew > rBest
                alpha = alpha / 2;
                rNew  = norm(F(qk - alpha * dx),2);
            end

        else
            alpha = 1;
        end

        % Update the nonlinear iterate in the (possibly bactracked) (quasi-) Newton direction
        qk     = qk - alpha * dx;
        rk     = F(qk)          ;
        Fqk    = rk             ;
        rBest  = norm(rk,2)     ;
        Q(:,1) = rk./rBest      ;


        % Perform convergence checks and iteration stepping
        NotDone             = norm(rk,2) > Tolerance    ;
        IterationsNonlinear = IterationsNonlinear + 1   ;
    end

end
