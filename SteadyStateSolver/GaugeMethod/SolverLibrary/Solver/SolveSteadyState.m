function [qSol,varargout] = SolveSteadyState(q)
    
    Requests = {'FluxSystemConnectivity','SourceSystemConnectivity','Ncv','dx'};
    [ConnectivityF,ConnectivityS,Ncv,dx] = Database('get',Requests{:});
    
    Requests = {'Flux','Source','FluxJacobian','SourceJacobian'};
    [Flux,Source,FluxJacobian,SourceJacobian] = Database('get',Requests{:});
    
    Requests = {'InitializationHook','AfterUpdateHook'};
    [InitializationHook,AfterUpdateHook] = Database('get',Requests{:});
    
    Gauge1 = [zeros(1,0*Ncv),1,zeros(1,3*Ncv-1)];
    Gauge2 = [zeros(1,1*Ncv),1,zeros(1,2*Ncv-1)];
    Gauge3 = [zeros(1,2*Ncv),1,zeros(1,1*Ncv-1)];
    
    if abs(det(ConnectivityF'*ConnectivityF)) > 100*eps()
        GaugedConnectivtyF = ConnectivityF;
        GaugedConnectivtyS = ConnectivityS;
        IsSingular = false;
    else
        %Gauge1;Gauge2;Gauge3;
        Gauged = [1,Ncv+1,2*Ncv+1];
        NotGauged = setdiff(1:Ncv*3,Gauged);
        GaugedConnectivtyF = ConnectivityF;
        GaugedConnectivtyS = ConnectivityS; 
        dx = dx(1);
        IsSingular = true;
    end
    
    qSol        = q     ;
    NotDone     = true  ;
    Tolerance   = 1E-10 ;
    IterMax     = 500   ;
    Iter        = 0     ;
    
    switch(IsSingular)
        case true
            
            InitializationHook(q);
            
            while NotDone
                %                 BeforeUpdateHook(q);
                
                Fluxes          = Flux(q)           ;
                Sources         = Source(q)         ;
                FluxJacobians   = FluxJacobian(q)   ;
                SourceJacobians = SourceJacobian(q) ;
                
                if Iter == 0
                    Scale = diag(1./Fluxes);
                end
                
                Residuals         = Scale*(GaugedConnectivtyF*Fluxes           -...
                                    GaugedConnectivtyS*Sources)          ;
                ResidualJacobians = Scale*(GaugedConnectivtyF*FluxJacobians    -...
                                    GaugedConnectivtyS*SourceJacobians)  ;
                
                [Q,R] = qr([Gauge1;Gauge2;Gauge3;ResidualJacobians])   ;
                dq    = - R\(Q'*[0;0;0;Residuals])      ;
                
                figure(1);
                semilogy(abs(dq));
                figure(2)
                plot(q(1:Ncv));
%                 RezAug   = [0;0;0;Residuals]                       ;
%                 RezJacAug = [Gauge1;Gauge2;Gauge3;ResidualJacobians];
%                 dq = (RezJacAug'*RezJacAug)\(RezJacAug'*RezAug);
%                 dq = -ResidualJacobians\Residuals;
                q  = q + dq;
                
                 ErrorLinf      = norm(dq,Inf)                      ;
                AboveTolerance = ErrorLinf > Tolerance              ;
                BelowIterMax   = Iter < IterMax                     ;
                NotDone        = AboveTolerance && BelowIterMax     ;
                
                Iter = Iter + 1;
                
                AfterUpdateHook(q);
                q(Ncv+1) = Database('get','rhouNew');

            end
            
        case false
    end
    
    
    
end