clc;
clear('all');


v0 = +0.529939;

% Control volume allocations
Mult = 1;
N.down = Mult*5  ;
N.low  = Mult*3  ;
N.heat = Mult*10 ;
N.chim = Mult*50 ;
N.hi   = Mult*10 ;
N.cool = Mult*10 ;

% Section lengths
L.down = 5                  ;
L.low  = 1                  ;
L.heat = 2                  ;
L.chim = L.down - L.heat    ;
L.hi   = L.low/2            ;
L.cool = L.low/2            ;

[q,Thermo,Info] = LoopSolverV2(v0,N,L);


% J_F = repmat({zeros(3)},length(q),1);
% J_S = repmat({zeros(3)},length(q),1);
%
% for k = 1:length(q)-1
%     State.T          = Thermo.T(k)      ;
%     State.P          = Thermo.P(k)      ;
%     State.f          = Thermo.f(k)      ;
%     State.Re         = Thermo.Re(k)     ;
%     State.f_Re       = Thermo.f_Re(k)   ;
%     State.mu         = Thermo.mu(k)     ;
%     State.PhaseCheck = true             ;
%
%     J_F{k} = CalculateFluxJacobian(q(:,k),State);
%     J_S{k} = CalculateSourceJacobian(q(:,k),State,Info.CV(k));
% end
%
%
% Mask1 = 1 : N.down;
% Mask2 = (N.down+1) : (N.down+N.low);
% Mask3 = (N.down+N.low+1) : (N.down+N.low+N.heat+N.chim);
% Mask4 = (N.down+N.low+N.heat+N.chim+1) : (N.down+N.low+N.heat+N.chim+N.hi+N.cool);
%
% J_F1Left  = J_F{Mask1( 1 )};
% J_F1Right = J_F{Mask1(end)};
% J_F2Left  = J_F{Mask2( 1 )};
% J_F2Right = J_F{Mask2(end)};
% J_F3Left  = J_F{Mask3( 1 )};
% J_F3Right = J_F{Mask3(end)};
% J_F4Left  = J_F{Mask4( 1 )};
% J_F4Right = J_F{Mask4(end)};
%
% J1Avg = (J_F1Left + J_F4Right)/2;
% J2Avg = (J_F2Left + J_F1Right)/2;
% J3Avg = (J_F3Left + J_F2Right)/2;
% J4Avg = (J_F4Left + J_F3Right)/2;
%
% L1 = L.down;
% L2 = L.low;
% L3 = L.heat + L.chim;
% L4 = N.hi + N.cool;
%
% J_S1 = CellMean(J_S(Mask1));
% J_S2 = CellMean(J_S(Mask2));
% J_S3 = CellMean(J_S(Mask3));
% J_S4 = CellMean(J_S(Mask4));
% J_SAvg = CellMean(J_S);
%
%
% A= blkdiag(J_S1-J2Avg/L1,J_S2-J3Avg/L2,J_S3-J4Avg/L3,J_S4-J1Avg/L4);
%
% A((1-1)*3+(1:3),(4-1)*3+(1:3)) = J1Avg/L1;
% A((2-1)*3+(1:3),(1-1)*3+(1:3)) = J1Avg/L1;
% A((3-1)*3+(1:3),(2-1)*3+(1:3)) = J1Avg/L1;
% A((4-1)*3+(1:3),(3-1)*3+(1:3)) = J1Avg/L1;
%
% D4 = eig(A);
% disp(' ');
% D1 = eig(J_SAvg);

%%
figure(1);
clf;
subplot(2,1,1)
xlabel('Node Number');
ylabel('Temperature [K]');
axis([1,89,372.9,373.75]);
grid('on');
box('on');
hold('on');
fill([9,9,19,19,9]   ,[370,380,380,370,370],'r','EdgeColor','r','FaceAlpha',0.3);
fill([79,79,89,89,79],[370,380,380,370,370],'b','EdgeColor','b','FaceAlpha',0.3);
fill([63,63,70,70,63],[370,380,380,370,370],'g','EdgeColor','g','FaceAlpha',0.3);
plot(Thermo.T,'k-','LineWidth',1.7);
hold('off');

subplot(2,1,2)
plot(Thermo.rho,'k-','LineWidth',1.7);
xlabel('Node Number');
ylabel('Density [kg/m^3]');
axis([1,89,200,1000]);
grid('on');
box('on');
hold('on');
fill([9,9,19,19,9]   ,[0,2E3,2E3,0,0],'r','EdgeColor','r','FaceAlpha',0.3);
fill([79,79,89,89,79],[0,2E3,2E3,0,0],'b','EdgeColor','b','FaceAlpha',0.3);
fill([63,63,70,70,63],[0,2E3,2E3,0,0],'g','EdgeColor','g','FaceAlpha',0.3);
hold('off');


figure(2);
clf;
hold('on');
plot(Thermo.rho,Thermo.i,'w');
rhoPlot = linspace(0.97*min(Thermo.rho),1.02*max(Thermo.rho),500);
Tisos   = linspace(min(Thermo.T  ),max(Thermo.T  ),20 );
for k = 1:length(Tisos)
    iPlot   = InternalEnergy(rhoPlot,Tisos(k));
    plot(rhoPlot,iPlot,'LineWidth',1.4,'Color',0.3*[1,1,1]);
end
plot(Thermo.rho,Thermo.i,'b','LineWidth',1.75);
axis([300,980,4.1875E5,4.215E5]);
box('on');
xlabel('Density [kg/m^3]');
ylabel('Internal Energy [J/kg]');
title('Overall Thermodynamic Path: Two-Phase Regime Emphasis');
MakePNG('ThermodynamicPath_2PhiRegime');

