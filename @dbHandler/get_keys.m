function keys = get_keys(obj, varargin)
keycell = obj.db.keys;
keys = {};
exclude = 'lorem ipsum dolor'; % certain not to be included
if ~isempty(varargin)
    pattern = varargin{1};
else
    keys = obj.db.keys';
end

if length(varargin) > 1
    exclude = varargin{2};
end

if ~isempty(varargin)
    for i=1:length(keycell)
        if contains(keycell{i}, pattern) && ~contains(keycell{i},exclude)
            keys = [keys; keycell{i}];
        end
    end
end
end