
Nconn = 100;

% Connectivities
SetOrdered = [1:3:3*Nconn,(1:3:3*Nconn)+1,(1:3:3*Nconn)+2];
LawOrdered = 1:3*Nconn;
Bidiag     = sparse(diag(-ones(Nconn,1),0)+diag(ones(Nconn-1,1),1));Bidiag =blkdiag(Bidiag,Bidiag,Bidiag);

LawBidiag  = Bidiag;
SetBidiag  = LawBidiag(:,SetOrdered);


% Jacobians
J = repmat({sparse([0,1,0;1,1,1;1,1,1])},Nconn,1);
JacobianSet = blkdiag(J{:});

J = sparse(1:Nconn,1:Nconn,1);
JacobianLaw = [0*J,J,0*J;J,J,J;J,J,J];


figure(1);
spy(JacobianSet);
figure(2);
spy(JacobianLaw);

figure(3);
spy(SetBidiag);
figure(4);
spy(LawBidiag);

figure(5);
spy(SetBidiag * JacobianSet);
figure(6);
spy(LawBidiag * JacobianLaw);

