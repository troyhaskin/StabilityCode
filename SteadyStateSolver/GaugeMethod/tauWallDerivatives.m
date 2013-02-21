function [tau_rho,tau_rhov,tau_rhoi] = tauWallDerivatives(q)
    
    % Database call
    Requests = {'Irho','Irhou','Irhoi','T','mu','PhaseCheck','Re','f','dfdRe','Lchar'};
    [Irho,Irhou,Irhoi,T,mu,PhaseCheck,Re,f,f_Re,Lchar] = Database('get',Requests{:});
    
    % Conserved variables
    rho  = q(Irho)  ;
    rhou = q(Irhou) ;
    rhoi = q(Irhoi) ;
    
    % Primitives
    Int   = rhoi ./ rho;
    rhou2 = rhou .* abs(rhou);
    
    % Thermodynamics
    [mu_rho,mu_i] = Viscosity_TotalDerivative(rho,T,Int,PhaseCheck);
    
    % Reynolds number derivatives
    Re_rho  = -Re./mu .* mu_rho         ;
    Re_rhou =  Lchar./mu .* sign(rhou)  ;
    Re_i    = -Re./mu .* mu_i           ;
    
    % Shear stress derivatives w.r.t. conserved variables
    tau_rho  = 2 * rhou2./(Lchar.*rho.^3) .*(f_Re.*(rho.^2 .* Re_rho - rhoi.*Re_i) - rho.*f);
    tau_rhov = 2 * (abs(rhou).*f + f.*rhou.*sign(rhou) + rhou2.*f_Re.*Re_rhou)./(rho.*Lchar);
    tau_rhoi = 2 * (1./Lchar) .* rhou2./(rho.^2) .* f_Re .* Re_i;
    
end


