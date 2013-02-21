function [dQ,dtMax] = GetTimeStepUpdates(Q,CFLLimit,Neqn,MaskL,MaskR,Ninter,dx,Flux,Source,InterfaceJacobian)
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
    
    dQ = Q*0                    ; % Allocate the output vector dQ
    qL = Q(Columnify(MaskL'))   ; % Pull the interface's left  states
    qR = Q(Columnify(MaskR'))   ; % Pull the interface's right states
    
    fL  = Flux(qL)      ; % Calculate the right flux values: size(fL) = [Nq,1]
    fR  = Flux(qR)      ; % Calculate the left  flux values: size(fR) = [Nq,1]
    psi = Source(Q) ; % Calculate the interface source terms
    
    
%     InterfaceFlux = @(q) fR(Mask) - fL(Mask) - dx(k)*psi(Mask);
    Ahats = InterfaceFluxJacobianExact(qL,qR); % Calculate the interface's intermediate flux Jacobian
    
    dtMax = Inf;
    MaskHandle = @(n) Neqn*(n-1) + 1;
    
    % ================================================================================ %
    %                       Iterate over all of the interfaces                         %
    % ================================================================================ %
    J = [1:Ninter,1:Ninter];
    for k = J
        
        Mask = MaskHandle(k) + [0,1,2];
        
        Ahat = Ahats(Mask,:);
        [R,D] = eig(Ahat);
        
        beta = R \ (fR(Mask) - fL(Mask) - dx(k)*psi(Mask));
        Speeds = diag(D);
        Zwaves = bsxfun(@times,R,beta');
        
        GoingLeft  = Speeds<0;
        GoingRight = Speeds>0;
        
%         fprintf('Interface %G Wave directions:\n',k);
%         Show([GoingLeft,GoingRight]');

        dQL = sum(Zwaves(:,GoingLeft) ,2);
        dQR = sum(Zwaves(:,GoingRight),2);
        
        dQ(MaskL(k,:)) = dQ(MaskL(k,:))+dQL./dx(k);%-dQR./dx(k);
        dQ(MaskR(k,:)) = dQ(MaskR(k,:))+dQR./dx(k);%-dQL./dx(k);
        
%         N = 1;subplot(3,1,N);plot(dQ(N:3:end),'bo','MarkerFaceColor','b');
%         N = 2;subplot(3,1,N);plot(dQ(N:3:end),'bo','MarkerFaceColor','b');
%         N = 3;subplot(3,1,N);plot(dQ(N:3:end),'bo','MarkerFaceColor','b');
        dtMax = min([dtMax;abs(CFLLimit.*dx(k)./(Speeds+eps()))/2]);
    end
    
end

function Addition = EmptySafeAdd(Operand1,Operand2)
    if not(isempty(Operand1)) && not(isempty(Operand2))
        Addition = Operand1 + Operand2;
        
    elseif isempty(Operand2) && not(isempty(Operand1))
        Addition = Operand1;
        
    elseif isempty(Operand1) && not(isempty(Operand2))
        Addition = Operand2;
    else
        Addition = Operand1 + Operand2;
    end
end