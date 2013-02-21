function dQ = GetSourceUpdates(Q,~,~,~,~,~,~,~,Source,~)
    %GetTimeStepUpdates
    %   FunctionName: GetTimeStepUpdates
    %   Purpose:    Calculate the updates to the conserved quantities in each control volume
    %
    %   Variables:
    %       Input:
    %           o Q
    %               Represents: the value of the conserved quantities in each control volume
    %               Type/Class: Numeric (double) array
    %               Size:       [Nq,1]
    %
    %           o MaskL
    %               Represents: Indices of all states to the "left" of every interface
    %               Type/Class: Numeric (double) array
    %               Size:       [Ninter,Neqn]
    %
    %           o MaskR
    %               Represents: Indices of all states to the "right" of every interface
    %               Type/Class: Numeric (double) array
    %               Size:       [Ninter,Neqn]
    %
    %           o Ninter
    %               Represents: number of interfaces
    %               Type/Class: Numeric (integer/whole double) scalar
    %               Size:       [1,1]
    %
    %           o dx
    %               Represents: spatial distance used for source intergration at every interface
    %               Type/Class: Numeric (double) array
    %               Size:       [Ninter,1]
    %
    %           o Flux
    %               Represents: handle to the flux function for the system
    %               Type/Class: function handle
    %               Returns:    Numeric (double) array of size [Nq,1]
    %
    %           o Source
    %               Represents: handle to the source function for the system
    %               Type/Class: function handle
    %               Returns:    Numeric (double) array of size [Ninter,1]
    %
    %           o InterfaceJacobian
    %               Represents: handle to the interface's intermediate flux Jacobian function
    %               Type/Class: function handle
    %               Returns:    Numeric (double) array of size [Nq,Neqn]
    %
    %       Output:
    %           o dQ
    %               Represents: updates to the conserved quantities
    %               Type/Class: Numeric (double) array
    %               Size:       [Nq,1]
    
    % ================================================================================ %
    %                                        Setup                                     %
    % ================================================================================ %
    
    dQ = Source(Q)              ; % Allocate the output vector dQ

end

