tic;
clc;
clear('all');

N = 10         ;
% A = wilkinson(N); 
A = rand(N,N) + diag(ones(N,1))   ;
x = rand(N,1)   ;
b = A*x         ;

Nrestart = 100;

tic;
[xMine,ResidualMine] = GMRES(A,b,[],N);
toc

tic;
xDirect = A\b;
toc

% semilogy(0:length(ResidualMine)-1,ResidualMine/ResidualMine(1),'r');


