clc;
clear('all');

Ncv  = 50       ;
Neqn = 3        ;
Nq   = Neqn*Ncv ;

% Jacobian-like formation
B = zeros(Nq,5);
d = zeros(1,5);

B(0*Ncv+1:end,1) = rand((Neqn-0)*Ncv,1);    
d(1) = 0;

B(1:2*Ncv,2) = rand(2*Ncv,1);    
d(2) = -1*Ncv;

B(1:1*Ncv,3) = rand(1*Ncv,1);    
d(3) = -2*Ncv;

B(2*Ncv+1:end,4) = 0.1*rand((Neqn-2)*Ncv,1);    
d(4) = +2*Ncv;

B(1*Ncv+1:end,5) = 100*rand((Neqn-1)*Ncv,1);    
d(5) = +1*Ncv;

JF = spdiags(B,d,Nq,Nq);

% Bidiagonal formation
Bi = spdiags([ones(Ncv,1),-[ones(Ncv-1,1);0]],[0,-1],Ncv,Ncv);
Bi = blkdiag(Bi,Bi,Bi);

figure(1);  spy(JF);
figure(2);  spy((Bi*JF)'*Bi*JF);
figure(3);  spy(Bi*JF);

condest(JF);
condest(Bi);
condest(Bi*JF);
condest((Bi*JF)'*Bi*JF);
