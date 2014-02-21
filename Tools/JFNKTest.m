clc;
% clear('all');
b = rand(2,1);
F = @(x) [x(1)*x(1) + x(1)*x(2) ; x(2)*x(1) + x(2)*x(2) ] + sin(x) - b  ;
J = @(x) [2*x(1) + x(2) + cos(x(1)) ,          x(1)             ;...
                   x(2)             , x(1) + 2*x(2) + cos(x(2)) ];


tic;
[xNewton,RezNewton,Iterations] = NewtonExact([10;50],F,J);
toc;
tic;
[xJFNK,IterationsNonlinear] = JFNKHouseholder([10;50],F,1E-3);
toc;

fprintf('Newton Solution (%G iterations):\n',Iterations);
disp(xNewton);
fprintf('JFNK Solution (%G iterations)\n',IterationsNonlinear);
disp(xJFNK);
