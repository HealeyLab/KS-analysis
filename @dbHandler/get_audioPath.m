function ap = get_audioPath(obj, workingDirectory)
    audioPaths = obj.audioPaths;
    
    % simplify audiopaths to just 'md_'
    for i=1:length(obj.audioPaths)
        path = audioPaths{i};
        parts = split(path);
        mdBlank = parts{end}; % md_
        
        if contains(workingDirectory, mdBlank)
            ap = path;
        end
    end
    
    if ~exist('ap', 'var')
        throw("no audiopaths matched")
    end
end