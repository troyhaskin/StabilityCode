tic;
clc;
clear('all');

N = 400         ;
A = wilkinson(N); % rand(N,N) + diag(ones(N,1))   ;
x = rand(N,1)   ;
b = A*x         ;

[xMine,Residual1] = GMRES(A,b,[],50);
% xMATs = gmres(A,b,50,1E-6,N);

semilogy(1:50,Residual1);
% semilogy(1:N,Residual1,1:N,linspace(9E-7,9E-7,N));
