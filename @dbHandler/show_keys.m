function keys = show_keys(obj, varargin)
    if ~isempty(varargin)
        keycell = obj.db.keys;
        pattern = varargin{1};
        
        keys = {};
        
        for i=1:length(keycell)
            if contains(keycell{i}, pattern)
                keys = [keys; keycell{i}];
            end
        end
    else
        keys = obj.db.keys';
    end
end