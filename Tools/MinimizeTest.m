

F = @(x) [ x(1) - x(2) ; x(1) + cos(x(2)) ]   ; % System
g = @(x) x(1)^2 + x(2)^2 - 1        ; % Constraint
L = @(x) norm(F(x),2) + x(3)*g(x)   ; % Langrangian: introduced x(3) as a mutiplier
R = @(x) [F(x);L(x)];

x0 = [0.5;0.5];

[x,RunStats] = JFNKHouseholder(x0,F,1e-8);