function [Q,R] = QRDecomposition(A)
    
    [M,N] = size(A) ;
    R     = A       ;
    Q     = eye(M)  ;
    e     = [1;zeros(M-1,1)];
    
    % Select number of iterations to perform
    if (M == N)         % Square
        Nmax = N - 1;
    elseif (N < M)      % Rectangular (more rows than cols)
        Nmax = N;
    elseif (N > M)      % Rectangular (more cols than rows)
        Nmax = M-1;
    end

    % Householder iteration
    for k = 1:Nmax

        % Grab k-th colum of R
        x = R(k:M,k);
        
        % Compute Householder vector for Householer Projection P
        u = -sign(x(1))*norm(x)*e(1:(M-k+1)) - x;
        u = u / norm(u);
        
        % Apply Householder projection I - 2*u*u' without forming P to R
        R = [ R(1:(k-1),:) ; R(k:M,:) - 2*u*(u'*R(k:M,:))];
        
        % Apply Householder projection I - 2*u*u' without forming P to Q
        Q = [ Q(1:(k-1),:) ; Q(k:M,:) - 2*u*(u'*Q(k:M,:))];
    end

    %   Let P = PNmax*...*Pk*...P2*P1
    %       o The for-loop performed P*A to produce R in an iterative fashion
    %       o Currently P*A = R
    %       o Since Q from the decomposition is unitary, P*A = Q'*A = R
    %       o Therefore, Q = P' or, using the variables of the function Q = Q'.
    Q = Q';
    R = triu(R);
end