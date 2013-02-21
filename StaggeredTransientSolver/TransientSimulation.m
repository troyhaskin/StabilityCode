function Qpoi = TransientSimulation(Solver,System,ControlVolumes,MomentumCells,Miscellaneous)
    
    % ====================================================================================== %
    %                           Add Transient Solver Library                                 %
    % ====================================================================================== %
    LibraryPaths = AddLibraryToPath();
    
    
    
    % ====================================================================================== %
    %                     Check User Input and Insert into Database                          %
    % ====================================================================================== %
    CheckAndRegisterInput('Solver'         , Solver         );
    CheckAndRegisterInput('System'         , System         );
    CheckAndRegisterInput('ControlVolumes' , ControlVolumes );
    CheckAndRegisterInput('MomentumCells'  , MomentumCells  );
    
    if (nargin == 5)
        CheckAndRegisterInput('Miscellaneous',Miscellaneous);
    end
    
    
    
    % ====================================================================================== %
    %                           Grab some useful constants                                   %
    % ====================================================================================== %
    tStart      = Solver.TimeStart      ; % Simulation start time
    tEnd        = Solver.TimeEnd        ; % Simulation end   time
    dtWrite     = Solver.dTimeSave      ; % Time stride for saving data
    dtMax       = Solver.dtMax          ; % Maximum time step to take
    NcvUnknowns = System.NcvUnknowns    ; % Number of unknown fields in the control volumes
    NmcUnknowns = System.NmcUnknowns    ; % Number of unknown fields in the control volumes 
    Ncv         = System.Ncv            ; % Number of control volumes
    Nmc         = System.Nmc            ; % Number of momentum cells
    Nq          = System.Nq             ; % Number of total conserved quantites in the system
    
    
    
    % ====================================================================================== %
    %                     Form the interface masks from the indices                          %
    % ====================================================================================== %
    FromVolumes = MomentumCells.FromControlVolumes(:);
    ToVolumes   = MomentumCells.ToControlVolumes(:);
    MCIndices   = (1:length(FromVolumes))';
    
    % Insert to database
    Database('insert','FromVolumes', FromVolumes    ,...
                      'ToVolumes'  , ToVolumes      ,...
                      'MCIndices'  , MCIndices      ); 
                      
    
    
    
    % ====================================================================================== %
    %                         Generate time point-of-interest vector                         %
    % ====================================================================================== %
    tPOI = (tStart : dtWrite : tEnd)'   ; % Time points where the Q states will be saved
    Npoi = length(tPOI)                 ; % Total number of time points to save.
    
    % Insert to database
    Database('insert','TimePointsOfInterest',tPOI);




    % ====================================================================================== %
    %                              Build Connectivity Arrays                                 %
    % ====================================================================================== %
    BuildConnectionArrays(FromVolumes,ToVolumes,MCIndices,Nmc,Ncv)




    % ====================================================================================== %
    %                                   Allocate Solution vector                             %
    % ====================================================================================== %
    
    Q    = zeros(Nq,Npoi);
    
    for k = 1:NcvUnknowns
        Mask      = (Ncv*(k-1)+1) : (Ncv*k);
        Q(Mask,1) = ControlVolumes.(['q',num2str(k)]);
    end
    for k = 1:NmcUnknowns
        Mask      = NcvUnknowns*Ncv + ((Nmc*(k-1)+1) : (Nmc*k));
        Q(Mask,1) = MomentumCells.(['q',num2str(k)]);
    end
    
    
    % ====================================================================================== %
    %                                        Time March                                      %
    % ====================================================================================== %
    Qpoi = TimeMarcher(Q,tPOI,dtMax);
    % Run time marcher
    % Post-process
    
    
    % ====================================================================================== %
    %                           Remove Transient Solver Library                              %
    % ====================================================================================== %
    RemoveLibraryFromPath(LibraryPaths);
    
end

function LibraryPaths = AddLibraryToPath()
    CurrentFolder = pwd();
    LibraryPaths  = genpath([CurrentFolder,'\TransientSolverLibrary']);
    addpath(LibraryPaths);
end

function [] = RemoveLibraryFromPath(LibraryPaths)
    rmpath(LibraryPaths);
end

