function [P_rhoTotal,P_iTotal] = Pressure_TotalDerivatives(rho,T,Int,PhaseCheck)
    
    if (nargin < 4) || isempty(PhaseCheck)
        PhaseCheck = true;
    end
    
    P_rho = Pressure_DensityN(rho,T,PhaseCheck);
    P_T   = Pressure_Temperature(rho,T,PhaseCheck);
    i_rho = InternalEnergy_DensityN(rho,T,PhaseCheck);
    i_T   = InternalEnergy_TemperatureN(rho,T,PhaseCheck);
    T_rho = -i_rho ./ i_T;
    T_i   = 1./i_T;
    
    P_rhoTotal = P_rho + P_T.*(T_rho - Int./rho .* T_i);
    P_iTotal = 1./rho .* P_T .* T_i;
    
end