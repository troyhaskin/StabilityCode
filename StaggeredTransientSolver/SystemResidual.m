function [Rez,RezJac] = SystemResidual(qNew,qOld,dt)
    
    persistent Irho Irhoi Irhou Aflow Vol Cmc dx g Lchar D f Jbase
    
    if isempty(Irho)
        RetrievalStruct = Database( 'get'                   ,... % GET request
            'Irho','Irhoi','Irhou'  ,... % Indices for extracting single quantites from q vectors
            'Irho_rhou','Jrho_rhou' ,... % Jacobian insertion indices
            'FromVolumes'           ,... % CV indexes that MC consider 'from' (i.e., left)
            'ToVolumes'             ,... % CV indexes that MC consider 'to'   (i.e., right)
            'FlowArea'              ,... % Flow area of Momentum Cells
            'Volume'                ,... % Volume of Control Volume
            'MomentumConnectivity'  ,... % Matrix to assist in Momentum summation
            'PathLength'            ,... % Momentum cell integration length
            'Gravity'               ,... % Orientation-dependent gravity vector
            'CharacteristicLength'  ,... % Length used for friction calculation
            'DiffusionCoefficient'  ,... % Difussion coefficient for water diffusion into atmosphere
            'FrictionFactor'        ,... % Friction factor
            'BaseJacobian'          ,... % Sparse matrix with all possible nonzero entries of the Jacobian
            'HeatAdded'             ,... % Heating function for control volumes
            'CurrentTime'           ,... % Current simulation time
            'MomentumConnections'   );   
        
        Irho  = RetrievalStruct.Irho                ;
        Irhoi = RetrievalStruct.Irhoi               ;
        Irhou = RetrievalStruct.Irhou               ;
        Irho_rhou   = RetrievalStruct.Irho_rhou     ;
        Jrho_rhou   = RetrievalStruct.Jrho_rhou     ;
        FromVolumes = RetrievalStruct.FromVolumes   ;
        ToVolumes   = RetrievalStruct.ToVolumes     ;
        Qadd        = RetrievalStruct.HeatAdded     ;
        Time        = RetrievalStruct.CurrentTime   ;
        Aflow = RetrievalStruct.FlowArea            ;
        Vol   = RetrievalStruct.Volume              ;
        Cmc   = RetrievalStruct.MomentumConnectivity;
        dx    = RetrievalStruct.PathLength          ;
        g     = RetrievalStruct.Gravity             ;
        Lchar = RetrievalStruct.CharacteristicLength;
        D     = RetrievalStruct.DiffusionCoefficient;
        f     = RetrievalStruct.FrictionFactor      ;
        Jbase = RetrievalStruct.BaseJacobian        ;
        CVMCIndices = RetrievalStruct.MomentumConnections;
    end
    
    % New conserved values
    rhoNew  = qNew(Irho);
    rhoiNew = qNew(Irhoi);
    rhouNew = qNew(Irhou);
    
    % Get the pressures
    P = GetPressures(rhoNew,rhoiNew);
    
    % Positivity masks
    PosMom    = (rhouNew >= 0);
    NegMom    = (rhouNew <  0);
    
    % Masks
    DonorMask =     PosMom  .* FromVolumes +     NegMom  .* ToVolumes;
    ZeroMask  = not(PosMom) .* FromVolumes + not(NegMom) .* ToVolumes;
    
    % Donor quantities
    rhoD      = rhoNew (DonorMask);
    rhoiD     = rhoiNew(DonorMask);
    pD        = P      (DonorMask);
    
    % Mass contributions to rho conservation
    MassIncrements = (Cmc * (rhouNew .* Aflow))./Vol - D.*rhoNew;
    
    % Energy contributions to rhoi conservation
    EnergyIncrements = (Cmc*(rhouNew./rhoD .* (rhoiD + pD) .* Aflow))./Vol + Qadd(Time);
    
    % Momentum contributions
    rhoAvg             = (rhoNew(ToVolumes)+rhoNew(FromVolumes))/2 ;
    MomentumIncrements = (P(FromVolumes) - P(ToVolumes))./dx                  +...
                            g .* rhoAvg                                       -...
                            2 * f()./Lchar .* rhouNew.*abs(rhouNew) ./ rhoAvg  ;
    
    % ======================================================================== %
    %                               Form Residual                              %
    % ======================================================================== %
    F   = [MassIncrements;EnergyIncrements;MomentumIncrements];
    Rez = (qNew - qOld) - F*dt;
    
    
    
    % ======================================================================== %
    %                           Form Residual Jacobian                         %
    % ======================================================================== %
    dFdqNew = Jbase;
    Lrho_rho  = sub2ind(size(Jbase),Irho,Irho);
    Lrho_rhou = sub2ind(size(Jbase),Irho_rhou,Jrho_rhou);
    dFdqNew(Lrho_rho) = -D;
    dFdqNew(Lrho_rhou) = Jbase(Lrho_rhou) .* sign(rhouNew(CVMCIndices)).*Aflow(CVMCIndices)./Vol(Irho_rhou);
    
end

function P = GetPressures(rho,rhoi)

    i = rhoi ./ rho;
    [Tguess,PhaseCheck] = Database('get','Tguess','PhaseCheck');

    [T,PhaseCheck] = Temperature(rho,i,Tguess,PhaseCheck);
    Database('set','Tguess',T,'PhaseCheck',PhaseCheck);

    P = Pressure(rho,T,PhaseCheck);

end





