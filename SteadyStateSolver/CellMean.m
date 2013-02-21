function Mean = CellMean(Cell)
    
    N = length(Cell);
    Accumulator = 0;
    
    for k = 1 : N
        Accumulator = Accumulator + Cell{k};
    end
    
    Mean = Accumulator/N;
    
end