function F = HEMFlux(q)
    
    % Get pertinent data from the Database
    [P,Irho,Irhou,Irhoi] = Database('get','P','Irho','Irhou','Irhoi');
    
    % Assign conserved variables
    rho  = q(Irho );
    rhou = q(Irhou);
    rhoi = q(Irhoi);

    % Construct primitives
    u = rhou ./ rho;

    % Form flux vector
    F        = q * 0            ;
    F(Irho)  = rhou             ;
    F(Irhou) = u .*  rhou + P   ;
    F(Irhoi) = u .* (rhoi + P)  ;

end