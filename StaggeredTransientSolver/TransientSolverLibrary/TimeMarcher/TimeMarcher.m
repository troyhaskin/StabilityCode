function Q = TimeMarcher(Q,tPOI,dtMax)
    
%     TimeStepUpdateParameters = GetTimeStepUpdateParameters();
    [PreTimeMarch,PreTimeStep,PostTimeStep,PostTimeMarch] = GetUserTimeStepFunction();
    
    CurrentTime = tPOI(1);
    Npoi = length(tPOI);
    Qwork = Q(:,1);
    Database('insert','CurrentTime',CurrentTime);
    
    
    % ============================================================= %
    %                 Execute User pre-march function               %
    % ============================================================= %
    PreTimeMarch(Qwork);
    
    
    % ============================================================= %
    %                        Begin time march                       %
    % ============================================================= %
    for k = 2:Npoi;

        NotATimePointOfInterest = true;

        % Loop until a point-of-interest has been found ------------------------
        while NotATimePointOfInterest

            % Execute User pre-step function ------------------------
            PreTimeStep(Qwork);

            % Adjust timestep to coincide with points-of-interest ---
            NextTime = CurrentTime + dtMax;
            if(NextTime < tPOI(k))
                dt = dtMax;
            else
                dt = tPOI(k) - NextTimeMax;
                NotATimePointOfInterest = false;
            end

            % Solve System ------------------------------------------
            Qwork = SolveImplicitSystem(Qwork,dt);

            % Update current time -----------------------------------
            CurrentTime = CurrentTime + dt;
            Database('set','CurrentTime',CurrentTime);

            % Execute User post-step function -----------------------
            PostTimeStep(Qwork);

        end

        % Push the the Point of Interest Q into the array ----------------------
        Q(:,k) = Qwork;

        % Execute User post-march function -------------------------------------
        PostTimeMarch(Qwork);

    end
    % ============================================================= %
    %                         End time march                        %
    % ============================================================= %
    
    
end


function Parameters = GetTimeStepUpdateParameters()
    Requests = {'MaskLForInterfaceSweep','MaskRForInterfaceSweep','Ninter',...
                'dx','Flux','Source','InterfaceFluxJacobian'};
    [MaskL,MaskR,Ninter,dx,Flux,Source,InterfaceJacobian] = Database('get',Requests{:});
    Parameters = {MaskL,MaskR,Ninter,dx,Flux,Source,InterfaceJacobian};
end

function [PreTimeMarch,PreTimeStep,PostTimeStep,PostTimeMarch] = GetUserTimeStepFunction()
    UserFunctions = {'PreTimeMarchFunction','PreTimeStepFunction','PostTimeStepFunction','PostTimeMarchFunction'};
    [PreTimeMarch,PreTimeStep,PostTimeStep,PostTimeMarch] = Database('get',UserFunctions{:});
end

