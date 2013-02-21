clc;
clear('all');

% Control volume allocations
Ncv = 8;

% Section lengths
Ldown = 5                  ;
Llow  = 1                  ;
Lheat = 2                  ;
Lchim = Ldown - Lheat    ;
Lhi   = Llow/2            ;
Lcool = Llow/2            ;



% Total control volume count
% Ncv = Ndown + Nlow + Nheat + Nchim + Nhi + Ncool;

% Frictional lengths for each control volume
Lfric = [1,1,1,5,5,5,1,1,1];

% Various scalar parameters
r      = 0.0508     ; % Pipe radius
Dh     = 2*r        ; % Hydraulic diameter
dx     = Lfric      ; % Integration distance / control volume length
Aflow  = pi*r.^2    ; % Flow area
Volume = dx.*Aflow  ; % Control volume volume
Lchar  = Dh         ; % Characteristic length
EpsRel = 1E-5       ; % Relative roughness
fGuess = 0.01       ; % Friction factor guess

% Gravity multiplier for each control volume
g = [0,0,0,-9.81,-9.81,-9.81,0,0,0];

% Heat load
qAddTotal  = 10E3/3;
% qHeatPerCV = qAddTotal/Nheat ;
% qCoolPerCV = -qAddTotal/Ncool ;

% 	qAddFlux = qAdd./(2*pi*r*dx);
qAddVol  = qAddTotal./Volume;
qAdd = [0,0,0,qAddVol,qAddVol,qAddVol,0,0,0];

%   Define the initial state - for a closed loop, this acts as the "pin" state.
%       o P = 101325.002    [Pa]
%       o T = 300           [K]
rho0 = 9.96556935266E+2 ;   % [kg/m^3]
i0   = 1.12553224581E+5 ;   % [J/kg]
v0   = 0.1              ;   %[m/s]
rhou = rho0 * v0;
rhoi = rho0 * i0;


% ============================================================================ %
%               Make ControlVolume information structs                         %
% ============================================================================ %
Sizer = ones(Ncv,1);

ControlVolume.FrictionalLengths     = Lfric     ;
ControlVolume.dx                    = dx        ;
ControlVolume.Gravity               = g         ;
ControlVolume.VolumetricHeatAddtion = qAddVol   ;
ControlVolume.Dh                    = Dh     * Sizer;
ControlVolume.Lchar                 = Lchar  * Sizer;
ControlVolume.EpsRel                = EpsRel * Sizer;
ControlVolume.qAdd                  = qAddVol;

ControlVolume.q1 = rho0 * Sizer ;
ControlVolume.q2 = rhou * Sizer ;
ControlVolume.q3 = rhoi * Sizer ;


% ============================================================================ %
%                   Make System information structs                            %
% ============================================================================ %

System.ControlVolumesLeft   = [1,2,3,2,3,4,5,6,7];
System.ControlVolumesRight  = [2,3,4,5,6,7,6,7,8];
System.Flux           = @(q) HEMFlux(q);
System.Source         = @(q) HEMSource(q);
System.FluxJacobian   = @(q) HEMFluxJacobian(q);
System.SourceJacobian = @(q) HEMSourceJacobian  (q);
System.Ncv            = Ncv      ;
System.Neqn           = 3         ;

Solver.InitializationHook = @(q) InitializeDatabase(q);
Solver.AfterUpdateHook    = @(q) UpdateSystem(q);
Solver.IntegrationOrder = 4;

q = SolveSteadyStateClosedLoop(Solver,System,ControlVolume,[]);



