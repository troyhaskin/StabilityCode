clc;
clear('all');

u   = @(x,t) exp(cos(5*pi/2*x).*sin(t/2));
u_t = @(x,t) u(x,t).*(cos(t/2).*cos(5*pi/2*x)/2)        ;
u_x = @(x,t) u(x,t).*(-5*pi/2*sin(t/2).*sin(5*pi/2*x))  ;

% MMS source for the PDE
% %   u_t + (cos(u)^2)_x = S
S = @(x,t) u_t(x,t) - sin(2*u(x,t)) .* u_x(x,t);

%{
% System Work put on hold.
v   = @(x,t) sin(5*pi/2*(x-t)).*cos(pi*x/2);
v_t = @(x,t) -5*pi/2*cos(5*pi/2*(x-t)).*cos(pi*x/2)        ;
v_x = @(x,t) pi*(5/2*cos(pi/2*x).*cos(5/2*pi*(x-t))-pi/2*sin(pi/2*x).*sin(5/2*pi*(x-t)));

% MMS sources for the system
% %   u_t + (cos(v))_x = Su
% %   v_t + (sin(u))_x = Sv
% Su = @(x,t) u_t(x,t) - sin(v(x,t)).*v_x(x,t);
% Sv = @(x,t) v_t(x,t) + cos(u(x,t)).*u_x(x,t);
%}


Nt       = 100000                               ;
Ncells   = 100                                  ;
t        = linspace(0.5,10,Nt  )'                 ;
dts      = t(2:end) - t(1:end-1)                ;
xEdges   = linspace(-1,1,Ncells+1)'             ;
xCenters = (xEdges(2:end) + xEdges(1:end-1))/2  ;
dxs      = xEdges(2:end) - xEdges(1:end-1)      ;


ConnectivityFluxExplicit   = spdiags([ [-ones(Ncells-1,1)./dxs(1:end-1);0],...
                                      [0;+ones(Ncells-1,1)./dxs(1:end-1)]],[-1,0],Ncells,Ncells);
ConnectivitySourceExplicit = speye(Ncells,Ncells);
ConnectivitySourceExplicit(1,1) = 0;

ubar = xCenters*0;
for k = 1:Ncells
    ubar(k) = quad(@(xp) u(xp,t(1)),xEdges(k),xEdges(k+1)) / (xEdges(k+1) - xEdges(k));
end
q = ubar;
F = @(qp,xp,tp) cos(qp).^2     ;
S = @(qp,xp,tp) S(xp,tp)    ;


plot(xCenters,u(xCenters,t(1)),xCenters,q,'o');
axis([-1,+1,0,3]);
legend('u_{exact}','u_{num}');
for k = 1:Nt
    
    qL   = q(1:end-1);
    qR   = q(2: end );
    qMid = [0;(qL + qR)/2];
    cR   = [(qR < qL)*+1 + (qR <= qL)*-1 ; 0];
    cL   = [0 ; cR(1:end-1)];
    
    q = q + dts(k) *  S(q,xCenters,t(k)) - ...
            dts(k) * (F([q(1);qR],xCenters,t(k))  - F([q(Ncells);qL],xCenters,t(k)))/dxs(1);
    
    plot(xCenters,u(xCenters,t(k)),xCenters,q,'o');
    axis([-1,+1,0,3]);
    drawnow();
    
    if k == 50
        g = [];
    end
    
end

% plot(xEdges,u(xEdges,t(1)),xEdges,v(xEdges,t(1)));
% axis([-1,+1,-1,3]);
% legend('u(x,t)','v(x,t)');
% for k = 2:Nt
%     plot(xEdges,u(xEdges,t(k)),xEdges,v(xEdges,t(k)));
%     axis([-1,+1,-1,3]);
%     drawnow();
% end