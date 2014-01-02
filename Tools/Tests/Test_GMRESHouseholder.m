clc;
clear('all');
% 

% Data load
A = load('west0479');
A = A.west0479      ;
N = length(A)       ;

% Definitions
x         = ones(N,1)       ; % solution vector
x0        = zeros(N,1)      ; % guess vector
b         = full(sum(A,2))  ; % Right-hand side
r0        = b               ; % Initial residual (b - A*0)
Nrestarts = 20              ; % Number of Arnoldi vectors to save
Nmax      = 20              ; % Maximum number of iterations
Tolerance = 1E-12           ; % Relative tolerance cutoff
nu        = 0.90            ; % ASGMRES weight factor

% Non-Preconditioned solves
PreConditionerLeft  = @(v) v;
PreConditionerRight = @(v) v;
[~,rNone] = GMRESHouseholderCore(A,r0,x0,Nrestarts,Nmax,Tolerance,nu,PreConditionerLeft,PreConditionerRight);
[~,~,~,~,rMat] = gmres(A,b,[],Tolerance,Nmax);

% Plot
subplot(2,1,1);
semilogy(0:length(rNone)-1,rNone/rNone(1),'-bo',0:length(rMat)-1,rMat/rMat(1),'-r*',...
    'LineWidth',1.75,'MarkerFaceColor','b');
xlabel('Iteration count');
ylabel('Relative Residual [-]');
legend('GMRESHouseholderCore','Matlab');
grid('on');




% Incomplete-LU Preconditioned version
[L,U] = ilu(A,struct('type','ilutp','droptol',1e-6));
PreConditionerLeft  = @(v) mldivide(U,mldivide(L,v));
PreConditionerRight = @(v) v ;
[xNone,rNone] = GMRESHouseholderCore(A,r0,x0,Nrestarts,Nmax,Tolerance,nu,PreConditionerLeft,PreConditionerRight);
[xMat,~,~,~,rMat] = gmres(A,b,[],Tolerance,Nmax,L,U);

% Plot
subplot(2,1,2);
semilogy(0:length(rNone)-1,rNone/rNone(1),'-bo',0:length(rMat)-1,rMat/rMat(1),'-r*',...
    'LineWidth',1.75,'MarkerFaceColor','b');
xlabel('Iteration count');
ylabel('Relative Residual [-]');
legend('GMRESHouseholderCore','Matlab');
grid('on');
