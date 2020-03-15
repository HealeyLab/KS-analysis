function keys = get_keys(obj, varargin)
%%GET_KEYS
%   Returns all keys that match the pattern
keycell = obj.db.keys;
keys = {};
exclude = 'lorem ipsum dolor'; % certain not to be included
if isempty(varargin)
    keys = obj.db.keys';
else
    pattern = varargin{1};
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