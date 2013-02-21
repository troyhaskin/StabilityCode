clc;
clear('all');

% Control volume allocations
Mult = 1;
Ndown = Mult*50 ;
Nlow  = Mult*10 ;
Nheat = Mult*25 ;
Nchim = Mult*25 ;
Nhi   = Mult*5 ;
Ncool = Mult*5 ;

% Section lengths
Ldown = 5                  ;
Llow  = 1                  ;
Lheat = 2                  ;
Lchim = Ldown - Lheat    ;
Lhi   = Llow/2            ;
Lcool = Llow/2            ;



% Total control volume count
Ncv = Ndown + Nlow + Nheat + Nchim + Nhi + Ncool;

% Frictional lengths for each control volume
LfricDown = Ldown/Ndown * ones(Ndown,1) ;
LfricLow  = Llow /Nlow  * ones(Nlow ,1) ;
LfricHeat = Lheat/Nheat * ones(Nheat,1) ;
LfricChim = Lchim/Nchim * ones(Nchim,1) ;
LfricHigh = Lhi  /Nhi   * ones(Nhi  ,1) ;
LfricCool = Lcool/Ncool * ones(Ncool,1) ;
Lfric     = [LfricDown;LfricLow;LfricHeat;LfricChim;LfricHigh;LfricCool];

% Various scalar parameters
r      = 0.0508     ; % Pipe radius
Dh     = 2*r        ; % Hydraulic diameter
dx     = Lfric      ; % Integration distance / control volume length
Aflow  = pi*r.^2    ; % Flow area
Volume = dx.*Aflow  ; % Control volume volume
Lchar  = Dh         ; % Characteristic length
EpsRel = 1E-5       ; % Relative roughness
fGuess = 0.01       ; % Friction factor guess

%         % Spatial Coordinates
%     xLeftDown = 0 * ones(Ndown,1);
%     xLeftLow  = Dh + cumsum([0;LfricLow(1:end-1)]);
%     xLeftHeat = (xLeftLow(end)+LfricLow(end)) * ones(Nheat,1);
%     xLeftChim = (xLeftLow(end)+LfricLow(end)) * ones(Nchim,1);
%     xLeftHigh = xLeftChim(end) - cumsum(LfricHigh(1:end-1));
%     xLeftCool = xLeftHigh(end) - cumsum([LfricCool(1:end);Dh]);
%     xLeft     = [xLeftDown;xLeftLow;xLeftHeat;xLeftChim;xLeftHigh;xLeftCool];
%
%     xRightDown = xLeftDown + Dh;
%     xRightLow  = xLeftLow  + cumsum(LfricLow(1:end));
%     xRightHeat = xLeftHeat + Dh;
%     xRightChim = xLeftChim + Dh;
%     xRightHigh = xLeftHigh + cumsum(LfricHigh(1:end));
%     xRightCool = xLeftCool + cumsum([LfricCool(1:end);Dh]);
%     xRight     = [xRightDown;xRightLow;xRightHeat;xRightChim;xRightHigh;xRightCool];

% Gravity multiplier for each control volume
g = [  +9.81 * ones(Ndown,1) ;...
       +0.00 * ones(Nlow ,1) ;...
       -9.81 * ones(Nheat,1) ;...
       -9.81 * ones(Nchim,1) ;...
       +0.00 * ones(Nhi  ,1) ;...
       +0.00 * ones(Ncool,1) ];

% Heat load
qAddTotal  = 10E3;
qHeatPerCV = qAddTotal/Nheat ;
qCoolPerCV = -qAddTotal/Ncool ;
qAdd = [               zeros(Ndown,1) ;...
                       zeros(Nlow ,1) ;...
         qHeatPerCV *  ones(Nheat,1) ;...
                       zeros(Nchim,1) ;...
                       zeros(Nhi  ,1) ;...
         qCoolPerCV *  ones(Ncool,1) ];

% 	qAddFlux = qAdd./(2*pi*r*dx);
qAddVol  = qAdd./Volume;

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

System.ControlVolumesLeft   = (1:Ncv)' ;
System.ControlVolumesRight  = ([2:Ncv,1])';
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



