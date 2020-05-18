function add(obj)
%% ADD Adds a new entry to dbHandler
% detailed info goes here
workingDirectory = pwd;
stim_data_path = fullfile(workingDirectory, 'adc_data.mat');
adc_data = load(stim_data_path); % S.board_adc, S.adc_sr
%% open contingency data
if ~isempty(dir('contingencies.txt'))
    fid = fopen('contingencies.txt','r');
    filecontents = fscanf(fid,'%c');
    fclose(fid);
    filecontents = strsplit(filecontents, '\n');
    % is now an array of 1x2 arrays of filename and samples
    filecontents = cellfun(@(x) strsplit(x, '\t'),filecontents,'UniformOutput',false)';

    if length(filecontents{end}) == 1
        filecontents(end) = [];
    end

    len = sum(cell2mat(cellfun(@(z) str2double(z{2}), filecontents,'UniformOutput',false)));

    %% Make habitBool
    habitBool = false(len,1);
    % count is the total number of samples we have processed so far
    count = 1;
    
    keyword = 'habit';
    newKeyword = input('Input keyword for differentiating between contexts (default: ''habit'')\n','s');
    if ~isempty(newKeyword)
        keyword = newKeyword;
    end
    
    
    % for each row of the textfile,
    for i = 1:length(filecontents)
        % Get the length of the current file
        curr_size = str2double(filecontents{i}{2});

        % If it's habit, make this chunk true. If it's negative, it's already
        % false.
        if contains(filecontents{i}{1},keyword)
            % We are building the ladder as we climb it: count is the start of
            % where we want to make it true, curr_size is the end of where we
            % want to make it true.

            % matlab is 1-indexed, I have to subtract 1 here to avoid an off by
            % one error
            habitBool(count:count+curr_size-1) = true;
        end

        % update count
        count = count + curr_size;
    end
else
    
    % set it to the length of the habituation. ASSUMES THIS IS A
    % HABITUATION TRIAL, THAT YOU AREN'T ADDING SHORT FEM RETURNS AGAIN
    habitBool = true(size(adc_data.board_adc(1, :), 2), 1);
end
% use `~habitBool` to indicate song times
%% 
fclose('all');            
has_stim_markers = ~isempty(dir('*markers.txt'));

if has_stim_markers
    % timestamps are the onset
    habit_adc = adc_data;
    habit_adc.board_adc = habit_adc.board_adc(:,habitBool);
    
    [stim_timestamps, stim_identities, ~] = obj.extract_stim_timestamps(habit_adc);
end

% Get Waveforms
if isempty(dir('raw_filtered.dat'))
    % Prerec for extracting waveforms
    obj.filter_raw_data()
end
[spike_waveform_info, spike_timestamps, sp] = obj.getWaveFormsDriver();

% Contingency caching
example_s = obj.remove(obj.get_family_name(workingDirectory));
if ~isempty(example_s)
    depth = example_s.depth;
    context = example_s.context;
    for_wf_analysis = example_s.for_wf_analysis;
    hemisphere = example_s.hemisphere;
else
    % get depth
    depth = input('depth of probe (um)\n');

    % get context
    context = '';
    while ~(strcmp(context, 'habituation') || strcmp(context, 'random')...
            || strcmp(context, 'song') || strcmp(context, 'other') ||...
            strcmp(context, 'E2'))
        context = input('habituation, random, song, E2, or other?\n','s');
    end

    % get hemisphere
    if contains(workingDirectory, 'mdx') || contains(workingDirectory, 'mde') || contains(workingDirectory, 'mdk')
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
end
%% put info into struct
% TODO: add LFP and if it's before or after E2
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
    s.amplifier_sampling_rate = sp.sample_rate;
    s.spike_timestamps = spike_timestamps(sp.clu == spike_waveform_info(i).unit);
    s.adc_sampling_rate = adc_data.adc_sr;

    % Assign spikes for song condition to song condition spikes and same
    % for pback spikes
    s.songSpks = slice(obj, ~habitBool, s);
    s.pbackSpks = slice(obj, habitBool, s);

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



    % for microphone trace, only add if context is "song"
    % I'm adding the mic data to a new key so that it doesn't
    % take up too much space, so that it doesn't take too much
    % time loading and saving
    
    if strcmp(s.context, 'song')
    
    % this gets kept with the rest of the data, redundantly. Because
    % I'm lazy and it's basically size zero.
        obj.db(obj.get_family_name(key)) = adc_data.board_adc(3,~habitBool);
    end
    % when was this added?
    s.whenadded = datetime;
    obj.db(key) = s;
end
beep;
%             obj.save_db();
end

function output = slice(obj, arrayBool, s)
    %% SLICE
    % Should return perfectly sorted timestamps
    
    output = [];
    offset = 0; on = 0;
    
    while any(arrayBool)
        
        % get all indices above/at threshold
        inds = find(arrayBool==1);
        
        % first index above threshold
        on = inds(1);

        % get all indicies after first index above threshold below threshold (ie, the downtick)
        inds = find(arrayBool(on:end)==0);
        if isempty(inds)
            off = length(arrayBool);
        else
            off = inds(1);
        end
        
        range = [on off] / s.adc_sampling_rate * s.amplifier_sampling_rate;            
        
        % timestamps (tss) within that range (s.spiketiemstamps is in samples)
        % note: these sample are already zero-ed. Next step is bringing it
        % to the front, it needed.
        tss = obj.sliceTS(s.spike_timestamps,range);
        
        % - on: aligns with the beginning of the recording
        % + offset: pushes it past all previous recordings
        % In most cases, the while loop runs once, and offset is zero.
        tss = tss + offset;
        
        % add these timestamps to the output
        output = [output tss];
        
        % update the offset: add the length of the current region to the summed lengths of all the other "on" regions
        offset = offset+diff(range);
        
        % To move on to the next envelope, I need to remove the current
        % "trues" so that when I find the next above threshold point, it is
        % not the same one
        arrayBool(1:off) = false;
    end
    output = output'; % fits with rest of db
end
