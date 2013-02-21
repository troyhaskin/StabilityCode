function [] = BuildConnectionArrays(FromVolumes,ToVolumes,MCIndices,Nmc,Ncv)
    
    
    % ============================================================================ %
    %                               Master row vectors                             %
    % ============================================================================ %
    Irho  = 0          + (1:Ncv)';
    Irhoi = Irho (end) + (1:Ncv)';
    Irhou = Irhoi(end) + (1:Nmc)';
    
    
    
    % ============================================================================ %
    %       Control Volume-Momentum Cell (CVMC) connectivity for Density           %
    %             and Internal Energy Conservation vector generation               %
    % ============================================================================ %
    Cmom      = zeros(Ncv,Nmc)   ;
    CVsRhous  = zeros(2*Nmc,1)   ;
    Irho_rhou = zeros(2*Nmc,1)   ;
    Jrho_rhou = zeros(2*Nmc,1)   ;
    Start     = 1                ;
    
    for k = 1:Ncv
        CVMClinkFrom = (FromVolumes == k)           ; % CVMC From links
        CVMClinkTo   = (ToVolumes   == k)           ; % CVMC To   links
        CVMClink     = (CVMClinkFrom | CVMClinkTo)  ; % CVMC All  links
        Nlink        = nnz(CVMClink)                ; % Total link count
        
        End                  = Nlink - 1 + Start            ; % End stride for insertion
        CVsRhous (Start:End) = MCIndices(CVMClink)          ; % MC indices for connected CV
        Irho_rhou(Start:End) = k                            ; % Rows    for Density  in JacobianBase matrix
        Jrho_rhou(Start:End) = Irhou(MCIndices(CVMClink))   ; % Columns for Momentum in JacobianBase matrix
        Start                = End + 1                      ; % Increment the insertion pointer
        
        
        % Fill momentum connectivity matrix used for flux summation over MC contributions to CVs
        Cmom(k,CVMClinkFrom) = -1   ;
        Cmom(k,CVMClinkTo  ) = +1   ;
    end
    % Sparsify the momentum connectivity matrix (probably not needed, but you never know)
    Cmom = sparse(Cmom);
    
    % Mask out unused entries (allocation was worst case).
    CVsRhous   = CVsRhous (1:End)       ;
    Irho_rhou  = Irho_rhou(1:End)       ;
    Jrho_rhou  = Jrho_rhou(1:End)       ;
    
    % Momentum connection for internal energies is the same
    Irhoi_rhou = Irho_rhou + Irho(end)  ; % Rows    for Internal Energy in JacobianBase matrix
    Jrhoi_rhou = Jrho_rhou              ; % Columns for Momentum        in JacobianBase matrix
    
    
    
    
    % ============================================================================ %
    %       Control Volume-Control Volume (CVCV) connectivity for Density          %
    %             and Internal Energy Conservation vector generation               %
    % ============================================================================ %
    Irhoi_rhoi = zeros(2*Nmc,1)   ;
    Jrhoi_rhoi = zeros(2*Nmc,1)   ;
    Start      = 1                 ;
    
    for k = 1:Ncv
        CVCVlinkFrom = (ToVolumes   == k);
        CVCVlinkTo   = (FromVolumes == k);
        
        CVCVlinkIndices       = [k ; FromVolumes(CVCVlinkFrom) ; ToVolumes(CVCVlinkTo)] ;
        Nlink                 = length(CVCVlinkIndices)                                 ;
        
        End                   = Nlink - 1 + Start   ; % End stride for insertion
        Irhoi_rhoi(Start:End) = k                   ; % Rows    for Internal Energy in JacobianBase matrix
        Jrhoi_rhoi(Start:End) = CVCVlinkIndices     ; % Columns for Internal Energy in JacobianBase matrix
        Start                 = End + 1             ; % Increment the insertion pointer
        
    end
    
    % Mask out unused entries and shift into Master index schem
    Irhoi_rhoi = Irhoi_rhoi(1:end) + Irho(end);
    Jrhoi_rhoi = Jrhoi_rhoi(1:end) + Irho(end);
    
    % Mask out unused entries and shift into Master index schem
    Irhoi_rho = Irhoi_rhoi(1:end)               ;
    Jrhoi_rho = Jrhoi_rhoi(1:end) - Irho(end)   ;
    
    
    
    
    % ============================================================================ %
    %       Momentum Cell-Control Volume (MCCV) connectivity for Momentum          %
    %                        Conservation vector generation                        %
    % ============================================================================ %
    Irhou_rho = [  MCIndices  ; MCIndices] + Irhoi(end) ;
    Jrhou_rho = [ FromVolumes ; ToVolumes]              ;
    
    Irhou_rhoi = [  MCIndices  ; MCIndices] + Irhoi(end);
    Jrhou_rhoi = [ FromVolumes ; ToVolumes] + Irho(end) ;
    
    
    
    
    % ============================================================================ %
    %                             Form Jacobian Base                               %
    % ============================================================================ %
    IrhoM = [Irho ; Irho_rhou];
    JrhoM = [Irho ; Jrho_rhou];
    
    %         Self    Other CVs     Momentum From           Momentum To
    IrhoiM = [Irhoi ; Irhoi_rhou ; Irhoi_rho ];
    JrhoiM = [Irhoi ; Jrhoi_rhou ; Jrhoi_rho ];
    
    %         Self           rhoi CV contributions           rho CV contributions
    IrhouM = [Irhou ; Irhou_rhoi ;  Irhou_rho];
    JrhouM = [Irhou ; Jrhou_rhoi ;  Jrhou_rho];
    
    Imaster      = [IrhoM;IrhoiM;IrhouM]        ;
    Jmaster      = [JrhoM;JrhoiM;JrhouM]        ;
    JacobianBase = sparse(Imaster,Jmaster,1)    ; 
    
    SizeJacobian = size(JacobianBase);
    Lrho_rho     = sub2ind(SizeJacobian, Irho      , Irho      );
    Lrho_rhou    = sub2ind(SizeJacobian, Irho_rhou , Jrho_rhou );
    Lrhoi_rho    = sub2ind(SizeJacobian, Irhoi_rho , Jrhoi_rho );
    Lrhoi_rhoi   = sub2ind(SizeJacobian, Irhoi     , Irhoi     );
    Lrhoi_rhou   = sub2ind(SizeJacobian, Irhoi_rhou, Jrhoi_rhou);
    Lrhou_rho    = sub2ind(SizeJacobian, Irhou_rho , Jrhou_rho );
    Lrhou_rhoi   = sub2ind(SizeJacobian, Irhou_rhoi, Jrhou_rhoi);
    Lrhou_rhou   = sub2ind(SizeJacobian, Irhou     , Irhou     );
    
    Database('insert',...
                'BaseJacobian', JacobianBase,...
                'Irho' , Irho ,'Irhoi' , Irhoi , 'Irhou' , Irhou        ,... % Indices for q extraction
                'Lrho_rho'    , Lrho_rho    , 'Jrho_rho'  , Irho        ,... % rho  - rho Indices
                'Lrho_rhou'   , Lrho_rhou   , 'Jrho_rhou' , Jrho_rhou   ,... % rho  - rhou Indices
                'Lrhou_rho'   , Lrhou_rho   , 'Jrhou_rho' , Jrhou_rho   ,... % rhou - rho  Indices
                'Lrhou_rhoi'  , Lrhou_rhoi  , 'Jrhou_rhoi', Jrhou_rhoi  ,... % rhou - rhoi Indices
                'Lrhou_rhou'  , Lrhou_rhou  , 'Jrhou_rhou', Irhou       ,... % rhou - rhou Indices
                'Lrhoi_rho'   , Lrhoi_rho   , 'Jrhoi_rho' , Jrhoi_rho   ,... % rhoi - rho  Indices
                'Lrhoi_rhoi'  , Lrhoi_rhoi  , 'Jrhoi_rhoi', Irhoi       ,... % rhoi - rhoi Indices
                'Lrhoi_rhou'  , Lrhoi_rhou  , 'Jrhoi_rhou', Jrhoi_rhou  ,... % rhoi - rhou Indices
                'MomentumConnections' , CVsRhous,... % Returns momentums ordered by Irho_rhou
                'MomentumConnectivity', Cmom);       % Coefficient matrix for summation
    
end
