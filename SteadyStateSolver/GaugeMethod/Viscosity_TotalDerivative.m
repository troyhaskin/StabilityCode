function [mu_rhoTotal,mu_iTotal] = Viscosity_TotalDerivative(rho,T,i,PhaseCheck)
    
    if (nargin < 4) || isempty(PhaseCheck)
        PhaseCheck = true;
    end
    
    mu_T   = Viscosity_TemperatureN        (rho,T,PhaseCheck);
    mu_rho = Viscosity_DensityN            (rho,T,PhaseCheck);
    T_i    = 1./InternalEnergy_TemperatureN(rho,T,PhaseCheck);
    i_rho  = InternalEnergy_DensityN       (rho,T,PhaseCheck);
    T_rho  = -T_i .* i_rho;
    
    mu_rhoTotal = mu_rho + mu_T.*(T_rho - i./rho .* T_i);
    mu_iTotal   = 1./rho .* mu_T .* T_i;
        
end