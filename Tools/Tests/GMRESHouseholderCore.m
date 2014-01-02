function [x,Residuals] = GMRESHouseholderCore(A,r0,x0,Nrestarts,Nmax,Tolerance,nu,...
        PreConditionerLeft,PreConditionerRight)
    
    % ================================================================================== %
    %                                GMRESHouseholderCore                                %
    % ================================================================================== %
    %{
        GMRESHouseholderCore is an implementation of Adaptive Simpler GMRES (ASGMRES)
        using Householder reflections for orthogonalization.  Standard GMRES relies on
        updating a Hessenberg matrix (the product of A and the update basis directions)
        and then solving a least-squares problem via the Hessenberg's QR factorization
        (with updates done via Givens rotations); Simpler GMRES starts with a slightly
        different basis vector such that only the QR factorization is needed.  Adaptive
        Simpler GMRES adds the weighting factor nu to iteratively change the choice of
        basis vector for the next Simpler GMRES iteration depending on the residual's
        behavior.  This adaptivity results in ASGMRES being more stable and well-
        conditioned for linear solves.


        Ref:    Rozlo�n�k Jir�nek, Pavel, and Miroslav Rozlo�n�k. "Adaptive version of
                Simpler GMRES." Numerical Algorithms 53.1 (2010): 93-112.
    %}
    
    
    
    % ============================================================= %
    %                             Set-Up                            %
    % ============================================================= %
    
    % Length of columns and vectors
    N = length(r0);
    
    % 
    Niterate = (Nrestarts <= Nmax)*Nrestarts + (Nrestarts > Nmax)*Nmax;
    
    % Matrix allocation
    Z = zeros(N,Niterate)  ;   % Update's basis vectors
    H = zeros(N,Niterate)  ;   % Householder vectors for projections
    Q = zeros(N,Niterate)  ;   % Unitary matrix columns of QR decomposition
    R = zeros(N,Niterate)  ;   % Upper-triangular matrix for least squares problem
    
    % Vector allocation
    e         = [1 ; zeros(N-1,1)]  ; % Unit vector used in creation of Householder vectors
    alpha     = r0*0                ; % Vector of projected residuals
    Residuals = zeros(Nmax,1)       ; % All residuals from ASGMRES
    
    % Initial residuals
    rk      = PreConditionerLeft(r0)    ; % Iterate residual
    rkNorm  = norm(rk,2)                ; % Iterated initial residual
    
    % Convergence iteration setup
    NotDone      = rkNorm > Tolerance   ;
    Tolerance    = rkNorm * Tolerance   ;
    Iterations   = 0                    ;
    Residuals(1) = rkNorm               ;
    n            = 2                    ;
    x            = x0                   ;
    
    
    % ================================================================================= %
    %                              Convergence Iterations                               %
    % ================================================================================= %
    
    while NotDone
        
        % ---------------------------------------------------- %
        %                    First Arnoldi Step                %
        % ---------------------------------------------------- %
        
        % First basis vector
        Z(:,1) = rk / rkNorm ;
        
        % Compute A*z1 and store in R
        Z(:,1) = PreConditionerRight(Z(:,1));
        R(:,1) = PreConditionerLeft (A*Z(:,1))   ;
        
        % Compute Householder vector for R(:,1)
        h      = R(:,1);
        h      = -sign(h(1)) * norm(h,2) * e(1:N) - h;
        H(:,1) = h / norm(h,2);
        
        % Bring R into upper triangular form via projection
        R(:,1) = R(:,1) - 2 * H(:,1) * (H(:,1)'*R(:,1));
        
        % Get the first column of the unitary matrix
        Q(:,1) = e - 2 * H(:,1) * (H(:,1)'*e);
        
        % Residual update
        alpha(1) = Q(:,1)'*rk           ; % Projected residual
        rk       = rk - alpha(1)*Q(:,1) ; % True b - A*(x+dx) residual
        
        % Starting and current residual norms used for choice of next basis vector
        rkm1Norm = rkNorm       ;
        rkNorm   = norm(rk,2)   ;
        
        % Store new residual
        Residuals(n) = rkNorm   ;
        n            = n + 1    ;
        
        
        % ---------------------------------------------------- %
        %                Requested Arnoldi Steps               %
        % ---------------------------------------------------- %
        for k = 2:Niterate
            
            if rkNorm < nu*rkm1Norm
                Z(:,k) = rk/rkNorm;
            else
                Z(:,k) = Q(:,k-1);
            end
            
            % Compute and store A*zk in R
            w      = PreConditionerRight(Z(:,k));
            R(:,k) = PreConditionerLeft (A*w)   ;
            
            % Apply all previous projections to new the column
            for m = 1:k-1
                h        = H(1:N-m+1,m);
                R(m:N,k) = R(m:N,k) - 2*h*(h'*R(m:N,k));
            end
            
            % Get the next Householder vector
            h            = R(k:N,k)                                 ;
            h            = -sign(h(1)) * norm(h,2) * e(1:N-k+1) - h ;
            h            = h / norm(h)                              ;
            H(1:N-k+1,k) = h                                        ;
            
            %   Apply projection to R to bring it into upper triangular form;
            R(:,1:k) = [R(1:k-1,1:k) ; R(k:N,1:k) - 2 * h * (h'*R(k:N,1:k))];
            
            % Get the k-th column of the current unitary matrix
            Q(:,k) = [zeros(k-1,1) ; e(1:N-k+1) - 2*h*(h'*e(1:N-k+1))];
            for m = k-1:-1:1
                hm       = H(1:N-m+1,m);
                Q(m:N,k) =  Q(m:N,k) - 2*hm*(hm'*Q(m:N,k));
            end
            
            % Update residual
            alpha(k) = Q(:,k)'*rk               ;
            rk       = rk - alpha(k)*Q(:,k)     ;
            
            % Reassign previous residuals
            rkm1Norm = rkNorm;
            rkNorm   = norm(rk,2);
            
            % Store new residual
            Residuals(n) = rkNorm   ;
            n            = n + 1    ;
            
            % Solve least-squares problem
            if rkNorm < Tolerance
                break;
            end
            
        end
        
        % ---------------------------------------------------- %
        %                   Post-Arnoldi Steps                 %
        % ---------------------------------------------------- %
        
        % Update iteration count
        Iterations = Iterations + k;
        
        % while-loop condition update
        NotConverged = (rkNorm > Tolerance)         ;
        CanIterate   = (Iterations <  Nmax)         ;
        NotDone      = NotConverged && CanIterate   ;
        
        
        % Update iteration maximum 
        if NotDone && (Nmax < (Iterations + Niterate))
            Niterate = Nmax - Iterations;
        end
        
        
        % Update solution

        % Solve the least-squares problem
        yk = triu(R(1:k,1:k)) \ alpha(1:k);
    
        % Calculate actual shift of guess
        dx = PreConditionerRight(Z(:,1:k)*yk);
        x  = x + dx;
        
    end

    Residuals = Residuals(1:k+1);
end
