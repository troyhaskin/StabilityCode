function [qBest,Thermo] = SolveSteadyStateClosedLoop(q0,Info,Tolerance)
    
    % Allocations
    Sizer  = zeros(1,Info.Ncv+1);
    q      = zeros(length(q0),Info.Ncv+1);
    Thermo = struct('rho',Sizer,'i',Sizer,'P',Sizer,'T',Sizer,'mu',Sizer,'Re',Sizer,'f',Sizer,'f_Re',Sizer);
    
    % Initialization
    q(:,1)        = q0(:);
    Thermo.rho(1) = q0(1);          % [kg/m^3]
    Thermo.i(1)   = q0(3)/q0(1);    % [J/kg]
    Thermo.T(1)   = Temperature(Thermo.rho(1),Thermo.i(1));
    Thermo.P(1)   = Pressure(Thermo.rho(1),Thermo.T(1));
    Thermo.mu(1)  = Viscosity(Thermo.rho(1),Thermo.T(1));
    Thermo.Re(1)  = ReynoldsNumber(q(2,1),Thermo.mu(1),Info.CV(1).Lchar);
    Thermo.f(1)   = FrictionFactorConstant(0.0137);
    Thermo.f_Re(1) = 0;
    
    % 'Best' solution
    qBest     = q;
    BestError = Inf;
    
    % Loop Setup
    NotDone   = true;
    Iter      = 0;
    
    % Print starting message
    DiagMessage = 'Iter:%3G  PdiffRel:%+10.2E Vel:%9s to %+9.6f SweepTime:%+10.2E\n';
    fprintf(DiagMessage,0,0,'set',q(2,1)/q(1,1),0);
    
    while NotDone
        % Start timer
        Tstart = tic;
        
        % Sweep through the system
        [q,Thermo] = SolveSteadyStateSweeper(q,Thermo,Info);
       
        % Get pressure error; for a closed loop P(end) == P(1)
        DeltaP   = Thermo.P(end) - Thermo.P(1)  ;
        ErrorAbs = abs(DeltaP)                  ;
        ErrorRel = ErrorAbs/Thermo.P(1)         ;
        NotDone  = abs(ErrorRel) > Tolerance    ;
        
        % Store 'best' solution
        if ErrorAbs < BestError
            qBest     = q;
            BestError = ErrorAbs;
        end
        
        % Adjust velocities according the pressures
        if NotDone
            [q(2,1),Outcome] = AdjustVelocity(q(2,1),DeltaP);
        else
            Outcome = 'unchanged';
        end
        
        % Odds and ends
        Iter = Iter + 1;
        Tend = toc(Tstart);
        fprintf(DiagMessage,Iter,ErrorRel,Outcome,q(2,1)/q(1,1),Tend);
    end
    
end


% Velocity corrections ----------------------------------------
function [q2New,Outcome] = AdjustVelocity(q2,DeltaP)
    
    persistent q2Old DeltaPold
    
    switch(not(isempty(q2Old)))
        case(true)
            
            % Get new guess from a linear approximation of the 0
            q2New = q2 - DeltaP./(DeltaP-DeltaPold)*(q2-q2Old);
            
            % Update old values
            q2Old     = q2;
            DeltaPold = DeltaP;
            
            if q2New > q2
                Outcome = 'increased';
            else
                Outcome = 'decreased';
            end
            
        case(false)
            
            % Store these values for later iterations
            q2Old     = q2;
            DeltaPold = DeltaP;
            
            % Half of double the speed for the first iteration
            if DeltaP > 0
                q2New= 1.1*q2;
                Outcome = 'increased';
            else
                q2New   = 0.9*q2;
                Outcome = 'decreased';
            end
    end
end

