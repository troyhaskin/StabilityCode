clc;
clear('all');

Np = 5E2;
x  = linspace(-1,1,Np)';
dx = x(2) - x(1);
dt = 0.005;

I = sparse(1:Np,1:Np,1);

% Advection matrix
A = -diag(sparse(ones(Np-1,1)),-1) + diag(sparse(ones(Np,1)),0);
A(1,1)       = 0 ;
A(1,2)       = 0 ; % Specified
A(end,end-1) = -dx/dt ; % Symmetry
A(end,end)   = 0 ; % Symmetry

% Diffusion matrix
D = diag(sparse(ones(Np-1,1)),-1) - 2 * diag(sparse(ones(Np,1)),0) + diag(sparse(ones(Np-1,1)),+1);

D(1,1)       = 0 ;
D(1,2)       = 0 ; % Specified
D(end,end-1) = 2 ; % Symmetry

disp(dt/dx^2);

T0  = (1 + tanh(500*x))/2;
TL = @(t) 0.5 - 0.5*sin(pi*(t+1/2));%.*exp(-0.01*t);

TdiffOld = T0;
TadvcOld = T0;

Time = 0;

while true
    TdiffNew  = (I - dt./dx.^2 * D) \ [TL(Time+dt);TdiffOld(2:end)];
    TadvcNew  = (I + dt./dx    * A) \ [TL(Time+dt);TadvcOld(2:end-1);0];
    
    plot(x,TdiffNew,x,TadvcNew,'LineWidth',2);
    axis([-1,+1,0,1.0]);
    drawnow 
    TdiffOld = TdiffNew;
    TadvcOld = TadvcNew;
    Time     = Time + dt;
end




