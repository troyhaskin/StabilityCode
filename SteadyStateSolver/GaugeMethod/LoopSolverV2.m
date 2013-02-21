function [q,Thermo,Info] = LoopSolverV2(v0,N,L)
    
    % Control Volume counts
    Ndown = N.down  ;
    Nlow  = N.low   ;
    Nheat = N.heat  ;
    Nchim = N.chim  ;
    Nhi   = N.hi    ;
    Ncool = N.cool  ;
    
    % Section lengths
    Ldown = L.down  ;
    Llow  = L.low   ;
    Lheat = L.heat  ;
    Lchim = L.chim  ;
    Lhi   = L.hi    ;
    Lcool = L.cool  ;
    
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
    qAddTotal = 0E3;
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
    rho0 = 9.96556935266E+2;    % [kg/m^3]
    i0   = 1.12553224581E+5;    % [J/kg]
    
    if (v0 == 0) || isempty(v0)
        v0 = 0.1; %[m/s]
    end
    
    q0 = [rho0 ; rho0 * v0 ; rho0 * i0];
    
    % Make the Flux, Source, and Residual handles
    Info.Ncv = Ncv;
    
    LfricCell = num2cell(Lfric)   ;
    dxCell    = num2cell(dx)      ;
    qAddCell  = num2cell(qAddVol) ;
    gCell     = num2cell(g)       ;
    % Tw      = num2cell(Tw)      ;
    Info.CV = struct('Lfric',LfricCell,'dx',dxCell,'qAdd',qAddCell,'g',gCell,...
        'Dh',Dh,'Lchar',Lchar,'EpsRel',EpsRel,'Tguess',300,'fGuess',fGuess);
    Info.Space   = [0;cumsum(dx)];
    
    
    TsimStart = tic;
    fprintf('************ Simulation  Started   ************\n');
    [q,Thermo] = SolveSteadyStateClosedLoop(q0,Info,1E-6);
    fprintf('************ Simulation Completed  ************\n');
    Info.RunTime = toc(TsimStart);
end


