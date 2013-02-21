function [] = InitializeThermodynamicProperties(q)

    [Irho,Irhoi] = Database('get','Irho','Irhoi');
    
    rho  = q(Irho)  ;
    rhoi = q(Irhoi) ;
    i    = rhoi ./ rho;
    [T,PhaseCheck] = Temperature(rho,i);
    
    Database('insert','PhaseCheck',                    PhaseCheck         ); % PhaseCheck boolean allocation
    Database('insert', 'Tguess'   ,                         T             ); % Temperature allocation
    Database('insert', 'P'        ,            Pressure(rho,T,PhaseCheck) ); % Pressure allocation
    Database('insert', 'mu'       ,           Viscosity(rho,T,PhaseCheck) ); % Viscosity allocation
    Database('insert', 'beta'     , VolumetricExpansion(rho,T,PhaseCheck) ); % Volumetric Expansion allocation
    Database('insert', 'k'        , ThermalConductivity(rho,T,PhaseCheck) ); % Thermal Conductivity allocation

end