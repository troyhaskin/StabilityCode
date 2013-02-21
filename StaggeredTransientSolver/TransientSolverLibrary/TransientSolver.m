function [q,varargout] = TransientSolver(SolverSetup,System,ControlVolumes,Interfaces,Miscellaneous)
    
    % Check user input
    TransientSolverCheckAndRegisterInput('SolverSetup'   ,SolverSetup   );
    TransientSolverCheckAndRegisterInput('System'        ,System        );
    TransientSolverCheckAndRegisterInput('ControlVolumes',ControlVolumes);
    TransientSolverCheckAndRegisterInput('Interfaces'    ,Interfaces    );
    
    if (nargin == 5)
        TransientSolverCheckAndRegisterInput('Miscellaneous',Miscellaneous);
    end
    
    % Run time marcher
    % Post-process
    q = 0;
end