function add(obj)
%% ADD Adds subject data to database
% 
    fclose('all');            
    workingDirectory = pwd;
    has_stim_markers = ~isempty(dir('*markers.txt'));
    if has_stim_markers
        % timestamps are the onset
        [stim_timestamps, stim_identities, adc_sr] = obj.extract_stim_timestamps(pwd);
    end

    % Get Waveforms
    if isempty(dir('raw_filtered.dat'))
        % Prerec for extracting waveforms
        obj.filter_raw_data()
    end
    [spike_waveform_info, spike_timestamps, sp] = obj.getWaveFormsDriver();

    % get depth
    depth = input('depth of probe (um)\n');
    
    % get context
    context = '';
    while ~(strcmp(context, 'habituation') || strcmp(context, 'random')...
            || strcmp(context, 'song') || strcmp(context, 'other'))
        context = input('habituation, random, song, or other?\n','s');
    end

    % get hemisphere
    if contains(workingDirectory, 'mdx')...
            || contains(workingDirectory, 'mde')...
            || contains(workingDirectory, 'mdk')
        hemisphere = 'R';
    elseif contains(workingDirectory, 'mdy') ...
            || contains(workingDirectory, 'mda')...
            || contains(workingDirectory, 'mdb')...
            || contains(workingDirectory, 'mdc')...
            || contains(workingDirectory, 'mdd')...
            || contains(workingDirectory, 'mdi')
        hemisphere = 'L';
    else
        error('This subject is not listed as having a hemisphere')
    end

    % use this recording for waveform analysis?
    yn = '';
    while ~strcmp(yn, 'y') && ~strcmp(yn, 'n')
        yn = input('Should this recording be used for broad/narrow waveform analysis? y/n\n', 's');
        if yn == 'y'
            for_wf_analysis = 1;
        elseif yn == 'n'
            for_wf_analysis = 0;
        end
    end

    %% put info into struct
    % TODO: add LFP and if it's before or after E2
    obj.remove(workingDirectory);
    for i = 1:length(spike_waveform_info)

        key = obj.keyhash(workingDirectory, spike_waveform_info(i).unit,...
            spike_waveform_info(i).channel,...
            strtrim(spike_waveform_info(i).goodness));
        % In case the number of units changed, we just wipe the
        % slate clean here

        % Now we (re)add everything
        s = struct;
        s.unit = spike_waveform_info(i).unit;
        s.channel = spike_waveform_info(i).channel;
        s.folder = workingDirectory; % encapsulates day and such
        s.goodness = spike_waveform_info(i).goodness; % 'MUA' or 'good'

        % need each unit's timestamps, not all of them
        s.spike_timestamps = spike_timestamps(sp.clu == spike_waveform_info(i).unit);
        s.amplifier_sampling_rate = sp.sample_rate;

        % actual waveform storage
        s.spike_waveforms = spike_waveform_info(i).waveform;
        s.for_wf_analysis = for_wf_analysis;

        % misc. info
        s.depth = depth;
        s.context = context;
        s.hemisphere = hemisphere;

        % adc stuff
        if has_stim_markers
            s.stim_timestamps = stim_timestamps;
            s.stim_identities = stim_identities;
        end

        
        stim_data_path = fullfile(workingDirectory, 'adc_data.mat');
        S = load(stim_data_path); % S.board_adc, S.adc_sr
        s.adc_sampling_rate = S.adc_sr;
        
        % for microphone trace, only add if context is "song"
        % I'm adding the mic data to a new key so that it doesn't
        % take up too much space, so that it doesn't take too much
        % time loading and saving
        if strcmp(s.context, 'song')
        % this gets kept with the rest of the data, redundantly. Because
        % I'm lazy and it's basically size zero.

            
            obj.db(obj.get_family_name(key)) = S.board_adc(3,:);
        end


        % when was this added?
        s.whenadded = datetime;
        obj.db(key) = s;
    end
%             obj.save_db();
end

