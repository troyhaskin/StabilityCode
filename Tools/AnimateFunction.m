function [] = AnimateFunction(F,xVector,ParameterVector,xLimits,yLimits)
    
    N = length(ParameterVector);
    
    for k = 1:N
        plot(xVector,F(ParameterVector(k)));
        axis([xLimits(:);yLimits(:)]');
        drawnow();
        pause(0.1);
    end
    
end 