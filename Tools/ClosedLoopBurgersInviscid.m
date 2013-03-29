clc;
clear('all');

k     = 10;
u0    = @(x) exp(-(k*x).^2);
u0bar = @(x,dx) sqrt(pi)*(erf(k*(x+dx/2))-erf(k*(x-dx/2)))./(2*k*dx);

x    = linspace(-1,1,1E3);
xAvg = linspace(-1,1,100);
xMid = (xAvg(2:end) + xAvg(1:end-1))/2;
xDel =  xAvg(2:end) - xAvg(1:end-1)   ;
plot(x,u0(x),xMid,u0bar(xMid,xDel),'ro');


Ncv      = 1000                             ;
xEdges   = linspace(-1,1,Ncv+1)'            ;
xDeltas  = xEdges(2:end) - xEdges(1:end-1)  ;
xCenters = xEdges(1:end-1) + xDeltas/2      ;

% Build the connectivity matrix
Ones = ones(1,Ncv-1);
I = [2:Ncv , 2:Ncv      , 1 , 1   ];
J = [2:Ncv , 1:(Ncv-1)  , 1 , Ncv ];
S = [Ones  , -Ones      , 1 , -1  ];
C = sparse(I,J,S);

D = zeros(Ncv,1);
D((xCenters>-0.6)&(xCenters<-0.4)) = +1;
D((xCenters<+0.6)&(xCenters>+0.4)) = -1;

plot(D)

NotDone = true;
q = 1/10*ones(Ncv,1);

C      = C(2:end,1:end);
Cprime = C              ;
Cprime(end+1,1) = 1;

while NotDone

    R  = [C*q - D(2:end) ; 0];
    
    subplot(2,1,1);
    plot(q);
    subplot(2,1,2);
    plot(R);
    
    dq = -Cprime\R;
    q  = q + dq;
    
    Show(norm(dq,1)/norm(q,1),'%+25.16E');
end


