function show_keys(obj, varargin)
    if ~isempty(varargin)
        keycell = obj.db.keys;
        pattern = varargin{1};
        for i=1:length(keycell)
            if contains(keycell{i}, pattern)
                disp(keycell{i})
            end
        end
    else
        disp(obj.db.keys')
    end
end