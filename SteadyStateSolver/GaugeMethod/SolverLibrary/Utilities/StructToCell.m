function Cell = StructToCell(Struct,AppendCell)
    
    FieldNames    = fieldnames(Struct)  ;
    Nfields       = length(FieldNames)  ;
    Cell          = cell(1,2*Nfields)   ;
    Cell(1:2:end) = FieldNames          ;
    
    for k = 1:Nfields
        Cell{2*k} = Struct.(FieldNames{k});
    end
    
    if (nargin == 2)
        Cell = [Cell,AppendCell];
    end
    
end