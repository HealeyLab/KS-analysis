function [cellArrOut] = sort_by_file_timestamp(~, fns)
%SORT_BY_FILE_TIMESTAMP Sorts cell array by the timestamp
%   when outputting from, say, dbHandler.get_keys(), it's sorted by 
%   filename, but chronological data is in the middle of the
%   file, so this takes that data in the middle of the file and sorts the
%   output according to it instead of the filename
cellArrOut=cell(length(fns),1);
tsCell = cellArrOut;

for i=1:length(fns)
    fns_split = split(fns{i}, '_');
    tsCell{i} = [fns_split{end-2:end-1}]; % now we have a timestamp
end

[~,i] = sort(tsCell);
cellArrOut = fns(i);
end

