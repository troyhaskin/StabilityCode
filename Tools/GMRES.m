function [x,Residual,LastResidualVec] = GMRES(A,b,x0,Nmax,r0)
    
    N = length(b);
    
    if (nargin < 4) || isempty(Nmax)
        Nmax = N;
    end
    
    if (nargin < 3) || isempty(x0)
        x0 = zeros(N,1);
    end
    
    if (nargin < 5) || isempty(r0)
        r0 = b - A*x0;
    end
    
    rNorm    = norm(r0,2)           ;
    H        = sparse([],[],[],Nmax,Nmax - 1,(Nmax+3)*Nmax/2+1);
    V        = zeros(N,Nmax)        ;
    Z        = zeros(Nmax,1)        ;
    I        = speye(Nmax,Nmax)     ;
    beta     = [rNorm ; Z(1:Nmax)]  ;
    V(:,1)   = r0./rNorm            ;
    Residual = Z                    ;


    % Initial GMRES step
    vk     = A*V(:,1)           ;
    H(1,1) = V(:,1)'*vk         ;
    vk     = vk - V(:,1)*H(1,1) ;
    H(2,1) = norm(vk,2)         ;
    V(:,2) = vk / H(2,1)        ;

    % Initialize QR matrices
    Hsq   = H(1,1)^2 + H(2,1)^2                         ;
    gk    = [rNorm;0]                                   ;
    Omega = sparse( [ H(1,1) ,  H(2,1)  ;...
                      H(2,1) , -H(1,1)] ) / sqrt(Hsq)   ;
    R     = sparse(sqrt(Hsq))                           ;
    
    for k = 2:Nmax
        
        % Form
        vk = A*V(:,k);
        
        % Form the k-th column of the Hessenberg
        H(1:k,k) = V(:,1:k)'*vk             ;
        vk       = vk - V(:,1:k)*H(1:k,k)   ;
        H(k+1,k) = norm(vk,2)               ;
        V(:,k+1) = vk/H(k+1,k)              ;
        
        % Calculate almost-triangular R
        Omega = [Omega , Z(1:k) ; Z(1:k)' , 1];
        R     = Omega * H(1:k+1,1:k);

        rho1 = R( k ,k);
        rho2 = R(k+1,k);
        rho  = sqrt(rho1^2 + rho2^2);
        rot1 = rho1 / rho;
        rot2 = rho2 / rho;

        G = blkdiag(I(1:k-1,1:k-1),[rot1,rot2;-rot2,rot1]);
        Omega  =  G    * Omega                  ;
        R      = Omega * H(1:k+1,1:k)           ;
        R      = triu(R)                        ;
        gk     = Omega * [rNorm;Z(1:k)]         ;
        yk     = R \ gk                         ;
        dx     = [V(1:k,1:k)*yk ; zeros(N-k,1)] ;
        rk     = r0 - A*dx                      ;
        rkNorm = norm(rk,2)                     ;
        
        Residual(k) = rkNorm;
    end
    
    Residual = [norm(r0,2);Residual];
    x =  (x0 + dx);
    LastResidualVec = rk;
    
end
