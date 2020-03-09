% post hoc bug fix: didn't add stim info before. Hopefully will now!
function reduce_size_of_db(obj)
    db = obj.db;
    keys = obj.db.keys;
%     db_aux = containers.Map('KeyType','char','ValueType', 'any'); % contains info for extra data, like microphone
%     aux_path = 'C:\Users\danpo\Documents\db_aux.mat';
    for i = 1:length(keys)
        s = db(keys{i});
        if isfield(s, 'microphone')
            % delete take the mic data
                            

            mic = s.microphone;
            s = rmfield(s, 'microphone');
            
            
            % move mic data to new key
            
            db(keys{i})  = s;
            db(familyname(keys{i})) = mic;
        end
    end
    % save manually, and first, make a backup copy
    obj.db = db;
    %     save(obj.dbPath, 'db','-v7.3');
%     save(aux_path, 'db_aux', '-v7.3');
    function fname = familyname(k)
        fname = strsplit(k, '&');
        fname = fname{1};
    end
        
        
end