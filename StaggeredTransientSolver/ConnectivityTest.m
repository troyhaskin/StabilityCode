function [] = BuildConnectionArrays(FromVolumes,ToVolumes,MCIndices,Nmc,Ncv)

% ============================================================================ %
%                     Control Volume and Momentum Cell counts                  %
% ============================================================================ %
% cvFrom = [1,2,3,4,5,6,3]'           ;
% cvTo   = [2,3,4,5,6,1,6]'           ;
% Ncv    = max([cvFrom(:);cvTo(:)])   ;
% Nmom   = length(cvFrom)             ;
% Mom    = (1:Nmom)'                  ;


% ============================================================================ %
%                               Master row vectors                             %
% ============================================================================ %
Irho  = 0          + (1:Ncv)';
Irhoi = Irho (end) + (1:Ncv)';
Irhou = Irhoi(end) + (1:Nmom)';
I     = [Irho ; Irhoi ; Irhou];



% ============================================================================ %
%       Control Volume-Momentum Cell (CVMC) connectivity for Density           %
%             and Internal Energy Conservation vector generation               %
% ============================================================================ %
Cmom      = zeros(Ncv,Nmom)   ;
CVsRhous  = zeros(2*Nmom,1)   ;
Irho_rhou = zeros(2*Nmom,1)   ;
Jrho_rhou = zeros(2*Nmom,1)   ;
Start     = 1                 ;

for k = 1:Ncv
    CVMClinkFrom = (cvFrom == k)                ; % CVMC From links
    CVMClinkTo   = (cvTo   == k)                ; % CVMC To   links
    CVMClink     = (CVMClinkFrom | CVMClinkTo)  ; % CVMC All  links
    Nlink        = nnz(CVMClink)                ; % Total link count
    
    End                  = Nlink - 1 + Start    ; % End stride for insertion
    CVsRhous (Start:End) = Mom(CVMClink)        ; % MC indices for connected CV
    Irho_rhou(Start:End) = k                    ; % Rows    for Density  in JacobianBase matrix
    Jrho_rhou(Start:End) = Irhou(Mom(CVMClink)) ; % Columns for Momentum in JacobianBase matrix
    Start                = End + 1              ; % Increment the insertion pointer
    
    
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
Irhoi_rhoi = zeros(2*Nmom,1)   ;
Jrhoi_rhoi = zeros(2*Nmom,1)   ;
Start      = 1                 ;

for k = 1:Ncv
    CVCVlinkFrom = (cvTo   == k);
    CVCVlinkTo   = (cvFrom == k);
    
    CVCVlinkIndices       = [k ; cvFrom(CVCVlinkFrom) ; cvTo(CVCVlinkTo)];
    Nlink                 = length(CVCVlinkIndices)     ;
    
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


IrhoM = [Irho ; Irho_rhou];
JrhoM = [Irho ; Jrho_rhou];

%         Self    Other CVs     Momentum From           Momentum To 
IrhoiM = [Irhoi ; Irhoi_rhou ; Irhoi_rho ];
JrhoiM = [Irhoi ; Jrhoi_rhou ; Jrhoi_rho ];

%         Self           rhoi CV contributions           rho CV contributions            
IrhouM = [Irhou ; [  Mom    ; Mom ] + Irhoi(end) ; [  Mom    ; Mom ] + Irhoi(end) ];
JrhouM = [Irhou ; [ cvFrom  ; cvTo] + Irho(end)  ; [ cvFrom  ; cvTo]              ];

Imaster = [IrhoM;IrhoiM;IrhouM];
Jmaster = [JrhoM;JrhoiM;JrhouM];
A = sparse(Imaster,Jmaster,eps()*eps());



end
