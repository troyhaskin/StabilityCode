tic;
clc;
clear('all');

N = 100         ;
Data = load('west0479');
A = Data.west0479;
x = rand(length(A),1)   ;
b = A*x         ;

Show(condest(A))

Nrestart = 20;

tic;
[xGMRES1,R1] = GMRES(A,b);
% [xGMRES1,R1] = GMRES(A,b,xGMRES1,N);
toc

xMat = gmres(A,b,N,1E-6,N);

figure(1);
spy(A);
figure(2);
semilogy(0:length(R1)-1,R1,'o');



norm(xGMRES1 - xMat,2)
norm(xGMRES1 - x,2)
norm(x - xMat,2)

