clc;
% clear('all');
% b = rand(2,1);
F = @(x) [x(1)*x(1) + x(1)*x(2) ; x(2)*x(1) + x(2)*x(2) ] + sin(x) - b  ;
J = @(x) [2*x(1) + x(2) + cos(x(1)) ,          x(1)             ;...
                   x(2)             , x(1) + 2*x(2) + cos(x(2)) ];


tic;
[xNewton,RezNewton,Iterations] = NewtonExact([0;0],F,J);
toc;
tic;
[xJFNK,IterationsNonlinear] = JFNKHouseholder([0;0],F,1E-9);
toc;

fprintf('Newton Solution (%G iterations):\n',Iterations);
disp(xNewton);
fprintf('JFNK Solution (%G iterations)\n',IterationsNonlinear);
disp(xJFNK);
