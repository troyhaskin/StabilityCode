tic;
clc;
clear('all');

N = 400         ;
A = wilkinson(N); % rand(N,N) + diag(ones(N,1))   ;
x = rand(N,1)   ;
b = A*x         ;

Nrestart = 100;

tic;
[xMine,ResidualMine] = GMRES(A,b,[],150);
toc
tic;
[xMAT,~,~,~,resvec]  = gmres(A,b,150,1.5E-15);
toc

% semilogy(1:length(resvec),resvec/resvec(1),'k',...
%          1:length(ResidualMine),ResidualMine/ResidualMine(1),'r--');
     
semilogy(1:length(resvec),resvec/resvec(1),'k',...
         1:length(ResidualMine),ResidualMine/ResidualMine(1),'r--');

% [xMine,Residual1] = GMRES(A,b,[],   Nrestart);
% [xMine,Residual2] = GMRES(A,b,xMine,Nrestart);
% [xMine,Residual3] = GMRES(A,b,xMine,Nrestart);
% [xMine,Residual4] = GMRES(A,b,xMine,Nrestart);
% semilogy(1:Nrestart,Residual1,1:Nrestart,Residual2,1:Nrestart,Residual3,1:Nrestart,Residual4);
