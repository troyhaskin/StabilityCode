function [] = TransientSimulationSetUp(varargin)

    RequiredInputsTop = {'Flux','Source','MaskLeft','MaskRight','TimeStart','TimeEnd','ControlVolumes'};
    RequiredInputsCV  = {'rho','rhov','rhoi','Orientation','CenterX','CenterY'};
%     OptionalFieldsTop = {'InterfaceFunction','dTimeSave'};

    %
    % Do some work to account for the Control Volumes struct argument
    ControlVolumeMask   = strcmp('ControlVolumes',varargin);
    ControlVolumeStruct = varargin{[false,ControlVolumeMask(1:end-1)]};
    ControlVolumeCell   = StructToCell(ControlVolumeStruct);
    
    % De-structified input (i.e., pure key-value list)
    Inputs = [varargin(notStructs),ControlVolumeCell];

    % Instantiate registry
    Registry();
    
    % Check for required, top-level requirements    
    for k = 1:length(RequiredInputsTop)
        
        MatchRequired = strcmp(RequiredInputsTop{k},Inputs);
        
        if not(any(MatchRequired))
            error('TransientSolver:TransientSimulationSetUp:MissingRequiredInput',...
                    'Required argument ''%s'' is not present.',RequiredInputsTop{k});
        else
            MatchRequiredInsertion = [false,MatchRequired(1:end-1)];
            Registry('insert',RequiredInputsTop{k},Inputs{MatchRequiredInsertion});
        end
    end
    
    
    % Check for required, CV-level requirements    
    for k = 1:length(RequiredInputsCV)
        if not(any(strcmp(RequiredInputsCV{k},Inputs)))
            error('TransientSolver:TransientSimulationSetUp:MissingRequiredCVInput',...
                    'Required control volume field ''%s'' is not present.',RequiredInputsCV{k});
        else
            Registry('insert',RequiredInputsCV{k},Inputs{k+1});
        end
    end
    
    
    % Check for optional arguments and assign defaults
    for k = 1:length(RequiredInputsCV)
        if not(any(strcmp(RequiredInputsCV{k},Inputs)))
            error('TransientSolver:TransientSimulationSetUp:MissingRequiredCVInput',...
                    'Required control volume field ''%s'' is not present.',RequiredInputsTop{k});
        end
    end
    
    
end



