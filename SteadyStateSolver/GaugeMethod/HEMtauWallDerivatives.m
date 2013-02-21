function [tau_rho,tau_rhov,tau_rhoi] = HEMtauWallDerivatives(q)
    
    % Database call
    Requests = {'Irho','Irhov','Irhoi','T','mu','PhaseCheck'};
    [Irho,Irhov,Irhoi,T,mu,PhaseCheck] = Database('get',Requests{:});
    
    % Conserved variables
    rho  = q(Irho);
    rhov = q(Irhov);
    rhoi = q(Irhoi);   
    
    % Primitives
    v   = rhov ./ rho;
    Int = rhoi ./ rho;
    
    % Geometry and flow data
    Lchar = Info.Lchar  ;
    f     = State.f     ;
    Re    = State.Re    ;
    f_Re  = State.f_Re  ;
    
    % Thermodynamics
    mu     = State.mu;
    T      = State.T;
    [mu_rho,mu_i] = Viscosity_TotalDerivative(rho,T,Int,false);
    
    % Reynolds number derivatives
    Re_rho  = -Re/mu * mu_rho;
    Re_rhov = Lchar/mu;
    Re_i    = -Re/mu * mu_i;
    
    % Shear stress derivatives w.r.t. conserved variables
    tau_rho  = 2 * (1/Lchar) * v^2/rho * (f_Re*(rho^2*Re_rho - rhoi*Re_i) - rho*f);
    tau_rhov = 2 * (1/Lchar) * v * (2*f + rhov * f_Re * Re_rhov);
    tau_rhoi = 2 * (1/Lchar) * v^2 * f_Re * Re_i;
    
end


