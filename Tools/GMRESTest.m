tic;
clc;
clear('all');

N = 100         ;
% A = wilkinson(N); 
A = rand(N,N) + diag(ones(N,1))   ;
x = rand(N,1)   ;
b = A*x         ;

Nrestart = 20;

tic;
[xGMRES1,R1] = GMRES(A,b,[],N);
toc

tic;
[xGMRES2,R2] = GMRES2(A,b,[],N);
toc

semilogy(0:length(R1)-1,R1,'o',0:length(R2)-1,R2);


