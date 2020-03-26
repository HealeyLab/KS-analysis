function [t] = sort_by(obj, keys, field)
%% SORT_BY Allows you to sort the inputted keys by a particular field
% input is in the form out dbHandler.get_keys()'s output, field is a string
% of a fieldname

    for i=1:length(keys)
        curr = obj.db(keys{i});
        if isfield(curr, field)
            currVal = curr.whenadded;
            newrow={keys{i} curr.(field)};
            if ~exist('tsCell', 'var')
                tsCell = cell(1,2);
                tsCell(1,:) = newrow;
            else
                tsCell = [tsCell; newrow];
            end
        end
    end
    keys = tsCell(:,1); % filtered
    [fieldOut,inds] = sort(cellfun(@(X) X(1), tsCell(:,2))); % sort by field
    cellArrOut = keys(inds); % sort the keys
    t = table(keys, fieldOut);
end