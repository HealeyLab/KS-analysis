function keysOut = filterKeys(obj, keys, field, pattern)
filter = true(length(keys),1);
for i = 1:length(keys)
    s = obj.db(keys{i});
    % if no good,
    if ~isstruct(s) || ~isfield(s,field) || ~compare(s,field,pattern)
        filter(i) = false;
    end
end
keysOut = keys(filter);
end
function boolOut = compare(s,field,pattern)
if ischar(pattern)
    boolOut = strcmp(s.(field), pattern);
elseif islogical(pattern)
    boolOut = s.(field) == pattern;
else
    error('type match error')
end
end