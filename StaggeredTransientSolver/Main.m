clc;
clear('all');


Neqn  = 3;
Nvert = 50;
Nhorz = 10;
Ncv   = 2*(Nvert + Nhorz);
Nq    = Ncv * Neqn;


% ============================================================================ %
%                                   Geometry                                   %
% ============================================================================ %

Hvert  = 5    ; %[m]
Hhorz  = 1.2  ; %[m]
l      = 0.1  ; %[m]
w      = 0.1  ; %[m]
Aflow  = l*w  ; %[m]
EpsRel = 1E-5 ; %[-]

SizerVert  = ones(Nvert,1);
SizerHorz  = ones(Nhorz,1);

% Downcomer
DowncomerOriens = cellstr(repmat('Vertical',Nvert,1))   ;
DowncomerH      = Hvert/Nvert                           ;
DowncomerHs     = SizerVert * DowncomerH                ;
DowncomerTops   = flipud(cumsum(DowncomerHs))           ;
DowncomerBots   = DowncomerTops - DowncomerH            ;
DowncomerLeft   = 0.0                                   ;
DowncomerLefts  = SizerVert * DowncomerLeft             ;
DowncomerRights = DowncomerLefts + w                    ;

% Lower run
LowerRunOriens = cellstr(repmat('Horizontal',Nhorz,1))              ;
LowerRunH      = (Hhorz - 2*l) / Nhorz                              ;
LowerRunHs     = SizerHorz * LowerRunH                              ;
LowerRunTops   = SizerHorz * l                                      ;
LowerRunBots   = LowerRunTops - l                                   ;
LowerRunLefts  = cumsum([DowncomerRights(end);LowerRunHs(2:end)])   ;
LowerRunRights = LowerRunLefts + LowerRunH                          ;

% Riser
RiserOriens = DowncomerOriens       ;
RiserH      = Hvert/Nvert           ;
RiserHs     = SizerVert * RiserH    ;
RiserTops   = cumsum(RiserHs)       ;
RiserBots   = RiserTops - RiserH    ;
RiserLeft   = LowerRunRights(end)   ;
RiserLefts  = SizerVert * RiserLeft ;
RiserRights = RiserLefts + w        ;

% Upper run
UpperRunOriens = LowerRunOriens                                 ;
UpperRunH      = (Hhorz - 2*l) / Nhorz                          ;
UpperRunHs     = SizerHorz * UpperRunH                          ;
UpperRunTops   = SizerHorz * Hvert                              ;
UpperRunBots   = UpperRunTops - l                               ;
UpperRunRights = cumsum([RiserLefts(end);-UpperRunHs(2:end)])   ;
UpperRunLefts  = UpperRunRights - UpperRunH                     ; 

% Master Dimensions
Orients = [DowncomerOriens ; LowerRunOriens ; RiserOriens ; UpperRunOriens ];
Lengths = [DowncomerHs     ; LowerRunHs     ; RiserHs     ; UpperRunHs     ];
Tops    = [DowncomerTops   ; LowerRunTops   ; RiserTops   ; UpperRunTops   ];
Bottoms = [DowncomerBots   ; LowerRunBots   ; RiserBots   ; UpperRunBots   ];
Lefts   = [DowncomerLefts  ; LowerRunLefts  ; RiserLefts  ; UpperRunLefts  ];
Rights  = [DowncomerRights ; LowerRunRights ; RiserRights ; UpperRunRights ];
CenterX = (Lefts   + Rights)/2;
CenterY = (Bottoms + Tops  )/2;
Volumes = Lengths * Aflow;
% Volumes = [Volumes,Volumes,Volumes]';
% Volumes = Volumes(:);



% ============================================================================ %
%                                     State                                    %
% ============================================================================ %
%   Define the initial state - for a closed loop, this acts as the "pin" state.
%       o P = 101325.002    [Pa]
%       o T = 300           [K]
rho0  = 9.96556935266E+2;    % [kg/m^3]
drho0 = 2.43953926053E-2;
i0    = 1.12553224581E+5;    % [J/kg]
di0   = -4.8065884983E+1;
v0    = 0.01;

DowncomerRhos =                  linspace(rho0,rho0+drho0,Nvert)' ;
DowncomerMoms = DowncomerRhos .* linspace(v0,v0,Nvert)' ;
DowncomerEnrs = DowncomerRhos .* linspace(i0,i0+di0,Nvert)'     ;

LowerRunRhos =                 linspace(rho0+drho0,rho0+drho0,Nhorz)'   ;
LowerRunMoms = LowerRunRhos .* linspace(v0,v0,Nhorz)'   ;
LowerRunEnrs = LowerRunRhos .* linspace(i0+di0,i0+di0,Nhorz)'       ; 

RiserRhos =              linspace(rho0+drho0,rho0,Nvert)' ;
RiserMoms = RiserRhos .* linspace(v0,v0,Nvert)' ;
RiserEnrs = RiserRhos .* linspace(i0+di0,i0,Nvert)'     ;

UpperRunRhos =                 linspace(rho0,rho0,Nhorz)'   ;
UpperRunMoms = UpperRunRhos .* linspace(v0,v0,Nhorz)'   ;
UpperRunEnrs = UpperRunRhos .* linspace(i0,i0,Nhorz)'       ;

rhos  = [DowncomerRhos ; LowerRunRhos ; RiserRhos ; UpperRunRhos];
rhovs = [DowncomerMoms ; LowerRunMoms ; RiserMoms ; UpperRunMoms];
rhois = [DowncomerEnrs ; LowerRunEnrs ; RiserEnrs ; UpperRunEnrs];




% ============================================================================ %
%                           Control Volume Structs                             %
% ============================================================================ %
ControlVolumes.Left                 = Lefts     ;   %% Optional
ControlVolumes.Right                = Rights    ;   %% Optional
ControlVolumes.Top                  = Tops      ;   %% Optional
ControlVolumes.Bottom               = Bottoms   ;   %% Optional
ControlVolumes.CenterX              = CenterX   ;
ControlVolumes.CenterY              = CenterY   ;
ControlVolumes.Length               = Lengths   ;   %% Optional
ControlVolumes.Volume               = Volumes   ;
ControlVolumes.Orientation          = Orients   ;
ControlVolumes.q1                   = rhos      ;
ControlVolumes.q2                   = rhois     ;



% ============================================================================ % 
%                               Heat Function                                  %
% ============================================================================ % 
Hheat        = 2.5; %[m]
Hcool        = 0.5; %[m]

MaskHeated    = (Bottoms>0.1) & (Tops<(0.1+Hheat)) & ((1:Ncv)'>Nvert);
Nheated       = sum(MaskHeated);
HeatedVolumes = Volumes .* MaskHeated;
HeatedVolume  = sum(HeatedVolumes);
HeatWeights   = HeatedVolumes./HeatedVolume;

MaskCooled    = (Lefts>0.1) & (Rights<(0.1+Hcool)) & ((1:Ncv)'>(2*Nvert+Nhorz));
Ncooled       = sum(MaskCooled);
CooledVolumes = Volumes .* MaskCooled;
CooledVolume  = sum(CooledVolumes);
CoolWeights   = -CooledVolumes./CooledVolume;

EnergySourceWeights = MaskHeated.*HeatWeights + MaskCooled.*CoolWeights;

HeatLoad   = 10E3   ; % [W]
TimeDead   = 0      ; % [s]
TimeRamp   = 1000   ; % [s]
Qdotdot    = HeatLoad / TimeRamp;
Ramp       = @(t) (Qdotdot*(t-TimeDead) .* HeatWeights  +  ... % Heating portion
                   Qdotdot*(t-TimeDead) .* CoolWeights) .* ... % Cooling portion
                  (t>=TimeDead).*(t<=TimeRamp)                                     +  ... % Turn off after ramp
                  (HeatLoad .* (HeatWeights + CoolWeights)).*(t>(TimeRamp+TimeDead));   % Constant Value


% ============================================================================ % 
%                          Momentum Cell information                           %
% ============================================================================ % 

FromControlVolumes = [(1:Ncv)' ;     03 ];
ToControlVolumes   = [(2:Ncv)' ; 1 ; 15 ];

Nmc    = length(FromControlVolumes)                                 ;
Sizer   = ones(Nmc,1)                                               ;
dx      = CenterX(ToControlVolumes) - CenterX(FromControlVolumes)   ;
dy      = CenterY(ToControlVolumes) - CenterY(FromControlVolumes)   ;
dCenter = abs(dx) + abs(dy)                                         ;
rhoAvg  = (rhos(FromControlVolumes) + rhos(ToControlVolumes))/2     ;

MomentumCells.q1                   = rhoAvg .* 0.01             ;
MomentumCells.FromControlVolumes   = FromControlVolumes         ;
MomentumCells.ToControlVolumes     = ToControlVolumes           ;
MomentumCells.Gravity              = -9.81*(abs(dy) > 100*eps());
MomentumCells.PathLength           = dCenter                    ;
MomentumCells.FlowArea             = Aflow             * Sizer  ;
MomentumCells.CharacteristicLength = 4*Aflow/(2*(l+w)) * Sizer  ;




% ============================================================================ %
%                               System Information                             %
% ============================================================================ %
System.Ncv          = Ncv                       ;
System.Nmc          = Nmc                       ;
System.NcvUnknowns  = 2                         ;
System.NmcUnknowns  = 1                         ;
System.NqCV         = System.NcvUnknowns * Ncv  ;
System.NqMC         = System.NmcUnknowns * Nmc  ;
System.Nq           = System.NqCV + System.NqMC ;



% ============================================================================ %
%                               Solver Information                             %
% ============================================================================ %
Solver.TimeStart = 0        ;
Solver.TimeEnd   = 2*3600   ;
Solver.dTimeSave = 1.0      ;
Solver.dtMax     = 0.1      ;
%
% Optional fields::
% Solver.TimeOut  = [];
% Solver.CFLlimit = [];
% Solver.dtMin    = [];
% Solver.dtMax    = [];
% Solver.Nsaves   = [];



% ============================================================================ %
%                                   Misc. Items                                %
% ============================================================================ %
Misc.HeatAdded             = Ramp                   ;
Misc.RelativeRoughness     = EpsRel                 ;
Misc.DiffusionCoefficient  = [3E-5; zeros(Ncv-1,1)];
Misc.FrictionFactor        = @(Re,EpsRel) 0.0137    ;
Misc.PreTimeMarchFunction  = @(q) InitializeThermodynamicProperties(q); %% Optional
Misc.PreTimeStepFunction   = @(q) [] ; % UpdateThermodynamicProperties(q);     %% Optional
Misc.PostTimeStepFunction  = @(q) [] ;    %% Optional
Misc.PostTimeMarchFunction = @(q) [] ;    %% Optional


% ============================================================================ %
%                                  Solver Call                                 %
% ============================================================================ %
q = TransientSimulation(Solver,System,ControlVolumes,MomentumCells,Misc);

% Solution = TransientSolver

