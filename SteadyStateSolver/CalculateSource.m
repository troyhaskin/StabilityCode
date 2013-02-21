function [S,State] = CalculateSource(q,State,Info)
    
    rho     = q(1);
    Gravity = Info.g;
    qAdd    = Info.qAdd;
    [dPdxFric,State] = tauWall(q,State,Info);
    
    S = [                0          ;...
           Gravity * rho - dPdxFric ;...
                       qAdd         ];
    
end