function extract_portmanteau(obj)
% TO get to this step in analysis:
% add cell
% run dbHandler.waveform_analysis and mark all included cells
% run set_playback_syllables for each subject
% run dbHandler.waveform_connector on each pair of cells
% run dbHandler.extract_portmanteau to summarize the data!
latency_db = load('C:\Users\danpo\Documents\latency_db.mat');
latency_db = latency_db.latency_db;
field_cell = fields(latency_db);
Z = @(B, S) (mean(S) - mean(B)) / sqrt(var(S) + var(B) - covar(S,B));

for i = 1:length(field_cell)
    orig_field = field_cell{i};
    % convert back to the two fields
    field = orig_field(2:end);
    % First, add back in the Kilosort
    field = strrep(field, '__', '_Kilosort_');
    % next, split up to two keys
    keys = strsplit(field, '_good');
    key1 = keys{1};
    key2 = keys{2};
    % finally, replace the last two underscores with ampersands.
    inds = strfind(key1, '_');
    key1(inds(end-1:end)) = '&';
    inds = strfind(key2, '_');
    key2(inds(end-1:end)) = '&';
    % Now it'll match a key in the database. Get'er done.
    key1 = char(obj.get_keys(key1));
    subject = obj.get_subject_id(key1);
    key2 = char(obj.get_keys(key2));
    
    % Habituation. key1.
    % assume always mde
    % get the FR and Z for each syllable
    BOS_latencies = latency_db.(orig_field).BOS; % struct of syllables!
    all_syllables = fields(BOS_latencies);
    for type_of_syllable = 1:length(all_syllables)
        syllable_cell = BOS_latencies.(all_syllables{type_of_syllable}); %'A', 'B', ... 'intro'
        latency_cell = cell(length(syllable_cell),1);
        
        Z_cell = cell(length(syllable_cell),1);
        
        % get average latency AND the other thing        
        for k = 1:length(syllable_cell)
            curr_row = syllable_cell{k};
            curr_row = curr_row(curr_row > 0);
            if isempty(curr_row)
                Z_cell{k} = 0;
                latency_cell{k} = NaN;
            else
                Z_cell{k} = length(curr_row);
                latency_cell{k} = curr_row(1);
            end
        end
        % latencies
        BOS_latencies.latency_cell = latency_cell;
        BOS_latencies.Z_cell = Z_cell;
        
        % Z scores
        B = obj.get_baseline(key1, 'BOS.wav');
        syl_arr = obj.get_playback_syllables(subject).BOS;
        % TODO: SORT SYL_ARR
        syl_arr = syl_arr(contains({syl_arr.id},'A'));
        
        % say there are 300 rows. There are syllables A, B, C, in the
        % order A, A, B, C, B, C. That means you need to find the order of
        % the syllables, and for that order, plus all the multples of that
        % order, get the zero-truncated window (the spikes have been
        % zero-truncated), and get the FR. Jesus.
        S = Z_cell{k} / diff(asdf);
        % get FR, transform to z score
    end
    
end
end