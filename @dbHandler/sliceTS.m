function output = sliceTS(obj, sp_ts, range)
        %% sliceTS, takes the cell of arrays rasterArr and 
        % returns only the values within the given time range.
        % is unit agnostic. Just make sure things are in the same unit.
        % inputs:
        %       sp_ts: timestamps
        %       range: range
        
        
        output = ((intersect(...
            sp_ts(sp_ts > range(1)),...
            sp_ts(sp_ts < range(2)))...
            - range(1)))';
    end