clc;
clear('all');


%   Need to load data for Tp, T, E0, E, R, F.
%   But I'm gonna make stuff up:
Tp = 290 + [1,2,3]      ;
T  = 300 + [10,20,30]   ;
E0 = [1,4,7]            ;
E  = [8,3,8]            ;
R  = 1                  ;
F  = 1                  ;


% Function handles
ConP = @(m) m(1) - m(2)./Tp ; % Plugging concentration as a function of fit parameters m = [a,b];
Con  = @(m) m(1) - m(2)./T  ; %          concentration as a function of fit parameters m = [a,b];
Sol  = @(m) E0 - R.*T./(2*F) .* log( ConP(m) ./ Con(m) ) ; % Solubility Curve


% Solve the least-squares problem.
Residual = @(m) E - Sol(m)              ; % Residual form used for fitting
mGuess   = [5,10]                       ; % Initial guess at coefficients
m        = lsqnonlin(Residual,mGuess)   ; % Solve the systems
% m = lsqnonlin(Residual,mGuess,mLowBound,mHighBound); % Bounds m's values


% Plot the resulting curve
plot(Tp,Sol(m),'b-o');



