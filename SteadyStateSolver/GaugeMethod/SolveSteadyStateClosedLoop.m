function [qSol,DatabaseSol] = SolveSteadyStateClosedLoop(Solver,System,ControlVolumes,Miscellaneous)
    
    % ====================================================================================== %
    %                           Add Transient Solver Library                                 %
    % ====================================================================================== %
    LibraryPaths = AddLibraryToPath();
    
    
    
    % ====================================================================================== %
    %                     Check User Input and Insert into Database                          %
    % ====================================================================================== %
    CheckAndRegisterInput([],Solver        );
    CheckAndRegisterInput([],System        );
    CheckAndRegisterInput([],ControlVolumes);
    CheckAndRegisterInput([],Miscellaneous );
    
    
    
    % ====================================================================================== %
    %                           Grab some useful constants                                   %
    % ====================================================================================== %
    AllocateAndFormCalculationArrays();
    Neqn   = System.Neqn            ; % Number of equations in the system
    Ncv    = System.Ncv             ;
    Nq     = Neqn * System.Ncv      ; % Number of total conserved quantites in the system
    
    
    
%     % ====================================================================================== %
%     %                     Form the interface masks from the indices                          %
%     % ====================================================================================== %
%     MaskL      = bsxfun(@plus,0:(System.Neqn-1),Neqn*(CVIndicesL-1)+1);
%     MaskR      = bsxfun(@plus,0:(System.Neqn-1),Neqn*(CVIndicesR-1)+1);
%     
%     % Insert to database
%     Database('insert','MaskLForInterfaceSweep',MaskL            ,...
%                       'MaskRForInterfaceSweep',MaskR            ,...
%                       'MaskL'                 ,Columnify(MaskL'),...  % Column versions of the sweep-friendly ones
%                       'MaskR'                 ,Columnify(MaskR'),...  % Column versions of the sweep-friendly ones
%                       'Ninter'                ,length(CVIndicesL));
    
    
    % ====================================================================================== %
    %                                   Allocate Solution vector                             %
    % ====================================================================================== %
    q = zeros(Nq,1);
    Bottom = 1;
    Top    = Ncv;
    for k = 1:Neqn
        q(Bottom:Top) = ControlVolumes.(['q',num2str(k)]);
        Bottom = Top + 1;
        Top    = Bottom + Ncv - 1;
    end


    % ====================================================================================== %
    %                                   Solve the System                                     %
    % ====================================================================================== %
    [qSol,DatabaseSol] = SolveSteadyState(q);
    
    
    % ====================================================================================== %
    %                           Remove Transient Solver Library                              %
    % ====================================================================================== %
    RemoveLibraryFromPath(LibraryPaths);
    
end

function LibraryPaths = AddLibraryToPath()
    CurrentFolder = pwd();
    LibraryPaths  = genpath([CurrentFolder,'\SolverLibrary']);
    addpath(LibraryPaths);
end

function [] = RemoveLibraryFromPath(LibraryPaths)
    rmpath(LibraryPaths);
end




