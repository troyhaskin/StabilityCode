function varargout = Database(Intent,varargin)
    
    % ============================================================================ %
    %                    Persistence and Initialization                            %
    % ============================================================================ %
    % Declare Database as persistent.  This makes permanent storage for
    % the variable and is only cleared when the base workspace is flushed
    % or a function is modified.
    persistent Database Options
    
    
    % If Database is newly created, it will be an empty array and will be
    % initialized as an empty struct.
    if isempty(Database)
        Database = struct();
        
        Options = struct();
        Options.InsertionsCannotOverwriteExistentFields = true;
    end
    
    
    
    % ============================================================================ %
    %                    Database-Wide Information/Diagnostics                     %
    % ============================================================================ %
    
    % Simple (i.e., no memory allocation) initialization
    if (nargin == 0)
        varargout{1} = true;
        return;
    end
    
    % Non-entry intents
    if (nargin == 1) && IsNotEmpty(Intent) && strcmp('_',Intent(1))
        switch(strtrim(lower(Intent)))
            case ('_fieldnames')
                varargout{1} = fieldnames(Database);
            case ({'_all','_registry','_database'})
                varargout{1} = Database;
            case ({'_whos','_information'})
                varargout{1} = whos('Database');
        end
        
        return;
    end
    
    % ============================================================================ %
    %                              Input Error checks                              %
    % ============================================================================ %
    ThrowingID = 'TransientSolver:Database:';
    
    if not(ischar(Intent))
        error([ThrowingID,'NonStringIntent'],...
            'The required input Intent must be empty or a string.');
    elseif not(ischar(varargin{1}))
        error([ThrowingID,'NonStringField'],...
            'The first requested field is not a string.');
    end
    
    
    
    % ============================================================================ %
    %                         Intent-dependent Set-Up                              %
    % ============================================================================ %
    
    % Semantic renaming of the data passed to the Database().
    IntentData = varargin;
    Nargout    = nargout();
    
    % Boolean Mask indicating which strings in IntentData are fields in the registry
    Mask       = isfield(Database,IntentData);
    
    % Pull the field names being dealt with from IntentData
    FieldNames = IntentData(Mask);
    
    % IntentDataMap Map:
    %   A vector whose entries indicate the index of the field name associated with
    %   the IntentData at the same index.
    %
    %   Example:
    %       o   Let IndexOptionMap = [1,1,2,2,2,3,3,3].
    %       o   This map says:
    %               * The first-second entries in IntentData are associated with the first  entry in FieldNames.
    %               * The third-fifth  entries in IntentData are associated with the second entry in FieldNames.
    %               * The sixth-eighth entries in IntentData are associated with the third  entry in FieldNames.
    %       o   Currently, the additional data is assumed to be indices for specific retrieval/assignment.
    %           In the future, non-numeric data may be used as options/switches, so the variable is named with
    %           that purpose in mind.
    %
    IntentData    = IntentData(not(Mask));    % IntentData without field names
    IntentDataMap = cumsum(Mask);             % Map with    field names
    IntentDataMap = IntentDataMap(not(Mask)); % Map without field names
    
    
    
    
    if any(Mask)
        % ============================================================================ %
        %                              Intent Switch                                   %
        % ============================================================================ %
        switch(lower(strtrim(Intent)))
            
            case {'retrieve','get'}
                GetRequestedDataFromDatabase();
                
            case {'assign','set'}
                AssignPassedDataToDatabase();
                varargout{1} = true;
                
            case 'delete'
                DeleteRequestedDataFromDatabase();
                varargout{1} = true;
                
            otherwise
                
        end
        
    else
        
        switch(strtrim(lower(Intent)))
            case 'initialize'
                InitializeDatabaseWithPassedData();
                varargout{1} = true;
                
            case 'insert'
                InsertPassedDataIntoDatabase();
                varargout{1} = true;
            otherwise
                
                % ************************* %
                %      INVALID REQUEST      %
                % ************************* %
                varargout{1} = false;
        end
        
    end
    
    
    
    % ============================================================================ %
    %                                Retrieval                                     %
    % ============================================================================ %
    
    % ----------------------- Main Routine ----------------------- %
    function [] = GetRequestedDataFromDatabase()
        
        Nfields = length(FieldNames);
        Output  = cell(1,Nfields);
        
        for k = 1:Nfields
            
            % Get the field name
            FieldName = FieldNames{k};
            
            % Get the registry entry
            Entry = Database.(FieldName);
            
            % Get the indices for the retrieval
            PossibleIndexData = IntentData(IntentDataMap == k);
            RetrievalIndices  = GetLinearIndices(PossibleIndexData,size(Entry));
            
            % Assign to output with/without indices
            if isempty(RetrievalIndices)
                Output{k} = Entry;
            else
                Output{k} = Entry(RetrievalIndices);
            end
            
        end
        
        if(Nargout >= Nfields)
            varargout = Output;
        else
            FieldDataPairs = cell(1,2*Nfields);
            FieldDataPairs(1:2:end) = FieldNames;
            FieldDataPairs(2:2:end) = Output;
            varargout{1} = struct(FieldDataPairs{:});
        end
        
    end
    
    
    
    % ============================================================================ %
    %                               Assignment                                     %
    % ============================================================================ %
    function [] = AssignPassedDataToDatabase()
        
        for k = 1:length(FieldNames)
            
            % Get the field name for assignment
            FieldName = FieldNames{k};
            
            % Get the entries full size
            EntrySize = size(Database.(FieldName));
            
            % Assignment Data
            AssignmentData   = IntentData(IntentDataMap == k);
            AssignmentValue  = AssignmentData{1};
            AssignmentIndices = GetLinearIndices(AssignmentData(2:end),EntrySize);
            
            % Assign to registry with/without indices
            if isempty(AssignmentIndices)
                % Without array subscripts, directly assign data.
                Database.(FieldName) = AssignmentValue;
                
            else
                % Without subscripts given.
                Database.(FieldName)(AssignmentIndices) = AssignmentValue;
                
            end
        end
    end
    
    
    % ============================================================================ %
    %                              Initialization                                  %
    % ============================================================================ %
    
    
    % ============================================================================ %
    %                              Insertion                                       %
    % ============================================================================ %
    function [] = InsertPassedDataIntoDatabase()
        
        FieldNames = IntentData(1:2:end);
        FieldData  = IntentData(2:2:end);
        Ninsert    = length(FieldNames);
        
        if isfield(Database,FieldNames)
            if Options.InsertionsCannotOverwriteExistentFields
                error([ThrowingID,'InsertionOverwriteField'],...
                    ['Insertion would overwrite an already existent field. ',...
                    'Changing existent values is done through setting. ']);
            end
        end
        
        for k = 1 : Ninsert
            Database.(FieldNames{k}) = FieldData{k};
        end
        
    end
    
    
    
    
    
    % ============================================================================ %
    %                                  Deletion                                    %
    % ============================================================================ %
    function [] = DeleteRequestedDataFromDatabase()
        Database = rmfield(Database,FieldNames);
    end
    
    
    % ============================================================================ %
    %                                  Helpers                                     %
    % ============================================================================ %
    
    function Indices = GetLinearIndices(DataSource,EntrySize_)
        
        % Indices should be numeric in addition to the Mask
        Mask       = cellfun(@isnumeric,DataSource);
        Subscripts = DataSource(Mask);
        
        if IsNotEmpty(Subscripts)
            
            try
                Indices = sub2ind(EntrySize_,Subscripts{:}); % Converts the subscripts to linear indices
            catch Error
                error('TransientSolver:Database:ImproperArraySubscripts'    ,...
                    ['The subscripts provided are invalid. This may be '  ,...
                    'a subscript that accidently exceeds the array''s '  ,...
                    'dimensions. Assignment/Set actions cannot expand '  ,...
                    'arrays.  If this is the source of the problem, '    ,...
                    'increase the initialized size of the the property.' ]);
            end
            
        else
            Indices = [];
        end
        
    end
    
    
    
end