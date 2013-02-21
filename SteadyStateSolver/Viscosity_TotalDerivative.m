function [mu_rhoTotal,mu_iTotal] = Viscosity_TotalDerivative(rho,T,Int,PhaseCheck)
    
    if (nargin < 4) || isempty(PhaseCheck)
        PhaseCheck = true;
    end
    
    mu_T   = Viscosity_TemperatureN     (rho,T,PhaseCheck);
    mu_rho = Viscosity_DensityN         (rho,T,PhaseCheck);
    i_T    = InternalEnergy_TemperatureN(rho,T,PhaseCheck);
    i_rho  = InternalEnergy_DensityN    (rho,T,PhaseCheck);
    T_rho  = - i_rho ./ i_T;
    T_i    = 1./i_T;
    
    mu_rhoTotal = mu_rho + mu_T.*(T_rho - Int./rho .* T_i);
    mu_iTotal   = 1./rho .* mu_T .* T_i;
        
end