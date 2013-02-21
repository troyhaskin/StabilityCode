function varargout = GetFromRegistry(varargin)
    
% Message ID for error/warning throws
    MessageID = 'TransientSolver:GetFromRegistry:';

% ---------------------------------------------------------------------------- %
%                Query the registry for the information                        %
% ---------------------------------------------------------------------------- %
    if iscellstr(varargin) && IsNotEmpty(varargin)
        RetrievedValues = Registry('retrieve',varargin);
    else
        error([MessageID,'NonStringInput'],...
              'The list of arguments to GetFromRegistry() must be strings only.');
    end


    
    
% ---------------------------------------------------------------------------- %
%               Return the data in the best way possible                       %
% ---------------------------------------------------------------------------- %
    Nout = nargout;
    Nin  = length(varargin);

    if ((Nout == 1) && (Nin > 1))
        % Return all requested values in a struct to the single output variable
        varargout{1} = RetrievedValues;
        
    elseif(Nout == Nin)
        % Return all requested values to each output variable one-to-one
        for k = 1:Nin
            varargout{k} = RetrievedValues.(lower(strtrim(varargin{k})));
        end
        
    elseif(Nout > Nin)
        % Return all requested vales to the pertinent outputs one-to-one, leave extras blank, and throw a warning.
        varargout{1:Nin} = RetrievedValues{1:Nin};
        warning([MessageID,'TooManyOutputs'],...
                ['The number of outputs was greater than the number requested. ',...
                 'Extraneous outputs variables are, therefore, empty.']);
        
    else
        % Return all requested values in a struct to the first output variable, leave extras blank, and throw a warning.
        varargout = cell2struct(RetrievedValues,varargin,1);
        warning([MessageID,'TooFewOutputs'],...
                ['The number of outputs was fewer than the number requested. ',...
                 'All retreived data was assigned to a struct in the first output ',...
                 'with all other outputs empty.']);
    end
    
    
end