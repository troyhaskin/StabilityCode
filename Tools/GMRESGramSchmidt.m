function [x,Residual,LastResidualVec] = GMRESGramSchmidt(A,b,x0,Nmax,r0)
    
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
%     H        = sparse([],[],[],Nmax,Nmax - 1,(Nmax+3)*Nmax/2+1);
%     V        = zeros(N,Nmax)        ;
%     Z        = zeros(Nmax,1)        ;
%     I        = speye(Nmax,Nmax)     ;
%     beta     = [rNorm ; Z(1:Nmax)]  ;
%     V(:,1)   = r0./rNorm            ;
%     Residual = Z                    ;
% 
% 
%     % Initial GMRES step
%     vk     = A*V(:,1)           ;
%     H(1,1) = V(:,1)'*vk         ;
%     vk     = vk - V(:,1)*H(1,1) ;
%     H(2,1) = norm(vk,2)         ;
%     V(:,2) = vk / H(2,1)        ;
% 
%     % Initialize QR matrices
%     Hsq   = H(1,1)^2 + H(2,1)^2                         ;
%     gk    = [rNorm;0]                                   ;
%     Omega = sparse( [ H(1,1) ,  H(2,1)  ;...
%                       H(2,1) , -H(1,1)] ) / sqrt(Hsq)   ;
%     R     = sparse(sqrt(Hsq))                           ;

    R      = zeros(N,N)     ;
    V      = zeros(N,Nmax)  ;
    Z      = zeros(Nmax,1)  ;
    Ar0    = A*r0           ;
    R(1,1) = norm(Ar0,2)    ;
    V(:,1) = Ar0 / R(1,1)   ;

    x      = x0         ;
    delta0 = rNorm      ;
    delta  = 1          ;
    r      = r0/rNorm   ;
    ksi    = Z(1:N)     ;
    rho    = Z(1:N)     ;
    
    for k = 2:Nmax
        
        % Form
        vk = A*V(:,k-1);
        
        for m = 1:k-1
            rho(m) = V(:,m)'*vk;
            vk     = vk - rho(m)*V(:,m);
        end
        rho(k) = norm(vk,2);
        V(:,k) = vk / rho(k);
        
        R(1:k,k) = rho(1:k);
        ksi(k)   = r'*V(:,k);
        delta    = delta*sin(acos(ksi(k)/delta));
        
        yk = R(1:k,1:k) \ ksi(1:k);
        
        z = yk(1)*r;
        for m = 2:k
           z = z + (yk(m) + yk(1)*ksi(m-1))*V(:,m-1);
        end
        x = x + delta0*z;
        
        r = (r - ksi(k)*V(:,k))/delta   ;
        delta0 = delta * delta0         ;
        delta  = 1                      ;
%         % Form the k-th column of the Hessenberg
%         H(1:k,k) = V(:,1:k)'*vk             ;
%         vk       = vk - V(:,1:k)*H(1:k,k)   ;
%         H(k+1,k) = norm(vk,2)               ;
%         V(:,k+1) = vk/H(k+1,k)              ;
% 
%         % Calculate almost-triangular R
%         Omega = [Omega , Z(1:k) ; Z(1:k)' , 1];
%         R     = Omega * H(1:k+1,1:k);
% 
%         rho1 = R( k ,k);
%         rho2 = R(k+1,k);
%         rho  = sqrt(rho1^2 + rho2^2);
%         rot1 = rho1 / rho;
%         rot2 = rho2 / rho;
% 
%         % Form Givens rotation matrix
%         G = [I(1:k-1,1:k-1) , Z(1:k-1) , Z(1:k-1) ;...
%               Z(1:k-1)'     ,   rot1   ,    rot2  ;...
%               Z(1:k-1)'     ,  -rot2   ,    rot1  ];
%         Omega  =  G    * Omega                  ;
%         R      = Omega * H(1:k+1,1:k)           ;
%         R      = (abs(R) > 100*eps()).*R        ;
%         gk     = Omega * [rNorm;Z(1:k)]         ;
%         yk     = R \ gk                         ;
%         dx     = [V(1:k,1:k)*yk ; zeros(N-k,1)] ;
%         rk     = r0 - A*dx                      ;
%         rkNorm = norm(rk,2)                     ;
%         
%         Residual(k) = rkNorm;
    end
    
    Residual = [norm(r0,2);Residual];
    x =  (x0 + dx);
    LastResidualVec = rk;
    
end
