% post hoc bug fix: didn't add stim info before. Hopefully will
% now!
function add_stim_info_to_all(obj)
    db = obj.db;
    keys = obj.db.keys;

    for i = 1:length(keys)
        s = db(keys{i});
        % first boolean, I removed the ~, making it positive
        disp(s.folder)

        if isfield(s,'stim_timestamps') && ~isempty(dir(fullfile(s.folder, '*markers.txt')))
            [stim_timestamps, stim_identities, adc_sr] = obj.extract_stim_timestamps(s.folder);
            s.stim_timestamps = stim_timestamps;
            s.stim_identities = stim_identities;
            s.adc_sampling_rate = adc_sr;
            obj.db(keys{i}) = s;

        end
    end
    save(obj.dbPath, 'db','-v7.3');
end