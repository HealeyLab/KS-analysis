%% Database handler with addition methods
% prop will point to .mat file with Containers.Map object comprised of
% db = containers.Map('KeyType','char','ValueType', 'any');
% will add to it with M(k) = v;

%% How to get this started
% db = containers.Map('KeyType','char','ValueType', 'any');
% save('db.mat', 'db','-v7.3')
%% How to use this
% First, review units in template-gui
% Second, run export_best_channels.py
% Third, find stimulus identities and times
% Finally, put everything into a database :(

% %
% x1=cursor_info_1.Position(1);
% x2=cursor_info_2.Position(1);
% plot(board_adc(2,:))
% board_adc(2,x1:x2) = 0;
% plot(board_adc(2,:))
% save('adc_data', 'board_adc', 'adc_sr', '-v7.3')
%%

classdef dbHandler
    properties
        dbPath = 'C:\Users\danpo\Documents\db_backup.mat';%'E:\DJP thesis sorting\db.mat';
        audioPathMDX = 'C:\Users\danpo\Documents\MATLAB\ephysSuite\zf son mdx';
        audioPathMDY = 'C:\Users\danpo\Documents\MATLAB\ephysSuite\zf son mdy';
        db = containers.Map; % initialized as none, basically
        count = 0;
        wf_keys = {} % a list of keys to use for wf_analysis
    end
    methods 
        function obj = dbHandler()
            S = load(obj.dbPath); 
            obj.db = S.db;   
            keys = obj.db.keys;
            for i = 1:length(keys)
                s = obj.db(keys{i});
                if s.for_wf_analysis
                    obj.wf_keys{length(obj.wf_keys) + 1} = keys{i};
                end
            end
        end

        % Make a key for the Map to point to the unit info
        function key = keyhash(~, workingDirectory, unit, channel, goodness)
            label = split(workingDirectory, '\');
            label = label{end};
            % key structure is dir&unit#&channel#, can split on '&'
            key = [label '&' num2str(unit) '&' num2str(channel) '&' goodness];
        end
        function [wd, unit, channel, goodness] = dehash(~, key)
            % returns info contained in keys
            splitcell = strsplit(key, '&');
            wd = splitcell{1};
            unit = splitcell{2};
            channel = splitcell{3};
            goodness = splitcell{4};
        end
            
        
        %% Remove all entries with match to key

        function remove(obj, key_family)
            key_family = strsplit(key_family, '&');
            key_family = key_family{1};
 
            db = obj.db;
            keys = db.keys;
            for i = 1:length(keys)
                if iscell(keys)
                    whole_key = keys{i};
                    key_cell = strsplit(whole_key, '&');
                    key_str = key_cell{1};
                    if strcmp(key_str, key_family)
                        remove(obj.db, whole_key);
                    end
                end
            end
        end
        
        function [on, off, wav_files] = get_stim_on_off(obj,key)
            s = obj.db(key);
            if contains(s.folder, 'mdx')
                audioPath = obj.audioPathMDX;
            else
                audioPath = obj.audioPathMDY;
            end
            
            wav_files = dir([audioPath '\*.wav']);
            for i = 1:length(wav_files)
                [cur_wav, fs] = audioread(fullfile(wav_files(i).folder,...
                    wav_files(i).name));
                cur_wav_resampled = resample(cur_wav, s.amplifier_sampling_rate, fs);
                wav_files(i).data = cur_wav_resampled(:,1);             
            end
            on = nan(length(s.stim_timestamps),1); off = nan(length(s.stim_timestamps),1);
            for i = 1:length(s.stim_timestamps)
                curr_wav_path = split(s.stim_identities{1}{i},'\');
                curr_wav_name = curr_wav_path{end};
                len = length(wav_files(strcmp({wav_files.name},curr_wav_name)).data);
                on(i) = s.stim_timestamps(i); %
                off(i)= s.stim_timestamps(i) + len; % / 2;
                
            end
        end
        
        function scatter(obj)
            vec=[];
            db = obj.db;
            for i = 1:length(obj.wf_keys)
                key = obj.wf_keys{i};
                s = db(key);
                disp(i)
                if strcmp(db(key).goodness, 'good')
                    %% So, my damn databse stores the center of the audio, not the beginning, so you need to find the beginning
                    if isfield(db(key), 'stim_timestamps')
                        [on, ~, ~]=get_stim_on_off(obj, key);

                        total_spikes = 0;
                        total_time = 0;
                        for j = 2:length(db(key).stim_timestamps)-1 % ignore the first and last one because might be truncated.
                            num_sp = length(s.spike_timestamps(...
                                s.spike_timestamps>(on(2)-s.amplifier_sampling_rate * 500) &...
                                s.spike_timestamps<on(2)));
                            total_spikes = total_spikes + num_sp;
                            total_time = total_time + 0.5;
                        end
                    else
                        spike_ts = db(key).spike_timestamps;
                        total_spikes = length(spike_ts);
                        total_time = spike_ts(end) - spike_ts(1);
                    end
                    %% 
                    vec = [vec; s.p2p s.sym total_spikes/total_time];
                end
            end
            scatter(vec(:,1), vec(:,2), 'filled','m')
            xlabel('p2p')
            ylabel('symmetry')
            
            figure;
            scatter3(vec(:,1), vec(:,2), vec(:,3))
            xlabel('p2p')
            ylabel('symmetry')
            zlabel('fr')
            title('Matt what do you think?')
            
        end
        function show_keys(obj)
            disp(obj.db.keys')
        end
        
        
        function keycell = cross_correlograms(obj)
            % First, separate all keys into 
            keys = obj.db.keys;
            keycell = cell(length(keys), 4 ); % rows, columns
            for i = 1:length(keys)
                [folder, unit, channel, goodness] = obj.dehash(keys{i});
                keycell{i,1} = folder;
                keycell{i,2} = unit;
                keycell{i,3} = channel;
                keycell{i,4} = goodness;
            end
            
            uniquecell = unique(keycell(:,1));
            
            % find cells from same recording
            for i = 1:length(uniquecell)
                recording_indices = find(strcmp(uniquecell{i},keycell(:,1)));
                
                % find cells from same tetrode
                recording_cells = keycell(recording_indices,:); % cells from same recording
                % note: channels can be 1 thru 16
                % for each tetrode, see if there are shared cells
                for j = 1:4
                    high = j * 4;
                    low = high - 3;

                    tetrode_cells = recording_cells(...
                        cellfun(@str2num, recording_cells(:,3)) >= low &...
                        cellfun(@str2num, recording_cells(:,3)) <= high,:);
                    
                    if size(tetrode_cells,1) > 1
                        obj.cross_corr_helper(tetrode_cells)
                    end
                end
            end
        end
        
        function cross_corr_helper(obj,tet_cells)
            db = obj.db;
            figure;
            hold on;
            for i = 1:size(tet_cells,1)
                key_i = tet_cells(i,:);
                s_i = db(obj.keyhash(key_i{1}, key_i{2}, key_i{3}, key_i{4}));
                
                for j = i+1:size(tet_cells,1)
                    
                    key_j = tet_cells(j,:);
                    s_j = db(obj.keyhash(key_j{1}, key_j{2}, key_j{3}, key_j{4}));
                    
                    [tsOffsets, ~, ~] = crosscorrelogram(...
                        s_i.spike_timestamps / 1000, s_j.spike_timestamps / 1000, [-0.100 0.100]);
                    
                    subplot(length(tet_cells),length(tet_cells),...
                        (i-1) * length(tet_cells) + j)
                    
                    hist(tsOffsets, 100);
                    xlim([-.1, .1]);
                end
            end      
            hold off;
        end
        
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
        %% Make an entry for each unit in the recording
        % Run per recording
        function add(obj) % to be called as dbHandler.add
            fclose('all');            
            workingDirectory = pwd;
            has_stim_markers = ~isempty(dir('*markers.txt'));
            if has_stim_markers
                [stim_timestamps, stim_identities, adc_sr] = obj.extract_stim_timestamps(pwd);
            end
            [spike_waveform_info, spike_timestamps, sp] = obj.getWaveFormsDriver();
            
            % get depth
            depth = input('depth of probe (um)\n');
            % get context
            context = '';
            while ~(strcmp(context, 'habituation') || strcmp(context, 'random')...
                    || strcmp(context, 'other')...
                    || strcmp(context, 't1') || strcmp(context, 't2')...
                    || strcmp(context, 't2') || strcmp(context, 't4'))
                context = input('habituation, random, or t1/t2/t3/t4, or other?\n','s');
            end
            
            % get hemisphere
            if contains(workingDirectory, 'mdy')
                hemisphere = 'L';
            elseif contains(workingDirectory, 'mdx')
                hemisphere = 'R';
            else
                hemisphere = '';
                while ~(strcmp(hemisphere, 'mdy') || strcmp(hemisphere, 'mdx'))
                    hemisphere = input('which hemisphere? (L/R)\n');
                end 
            end
            
            % use this recording for waveform analysis?
            yn = '';
            while ~strcmp(yn, 'y') && ~strcmp(yn, 'n')
                yn = input('Should this recording be used for waveform analysis? y/n\n', 's');
                if yn == 'y'
                    for_wf_analysis = 1;
                elseif yn == 'n'
                    for_wf_analysis = 0;
                end
            end
            
            % for each element in spike_waveform_info, which has info for
            % all units, 'good' and 'MUA'
            obj.remove(workingDirectory);
            for i = 1:length(spike_waveform_info)
                
                key = obj.keyhash(workingDirectory, spike_waveform_info(i).unit,...
                    spike_waveform_info(i).channel, spike_waveform_info(i).goodness);
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
                % complicated ones
                s.depth = depth;
                s.context = context;
                s.hemisphere = hemisphere;
                % adc stuff
                if has_stim_markers
                    s.stim_timestamps = stim_timestamps;
                    s.stim_identities = stim_identities;
                    s.adc_sampling_rate = adc_sr;
                end
                % when was this added?
                s.whenadded = datetime;
                obj.db(key) = s;
            end
            
            db = obj.db;
            save(obj.dbPath, 'db','-v7.3');
        end
        
        function [spike_waveform_info, spike_timestamps, sp] = getWaveFormsDriver(~)
            close all;
            
            dataDir = pwd;
            sp = loadKSdir(dataDir);
            
            gwfparams.sr = sp.sample_rate;
            gwfparams.dataDir = dataDir;
            gwfparams.fileName = 'raw_filtered.dat'; % this is so it errors out if you haven't filtetred the data already
            gwfparams.dataType = sp.dtype;
            gwfparams.nCh = sp.n_channels_dat;
            gwfparams.wfWin = [-(0.001*gwfparams.sr) 0.002*gwfparams.sr];  % Number of samples before and after spiketime to include in waveform
            % FRAC = 0.10; % fraction of all spikes to sample from
            
            gwfparams.nWf = 2000; % floor(length(sp.st) * FRAC);
            
            spike_timestamps = sp.st * sp.sample_rate;
            gwfparams.spikeTimes = ceil(spike_timestamps);
            
            gwfparams.spikeClusters = sp.clu;

            wf = getWaveForms(gwfparams);
            %%
            input('when you run the export best channels script, click enter')
            %%
            cluster_quality = readtable('cluster_groups.csv');
            best_channels   = readtable('best_channels.csv');
            good_clusters = cluster_quality(strcmp(cluster_quality.group, 'good'),1);
            MUA_clusters  = cluster_quality(strcmp(cluster_quality.group, 'mua'),1);
            good_clusters = table2array(good_clusters);
            MUA_clusters  = table2array(MUA_clusters);
            spike_waveform_info = []; % empty struct that will increase in size 
            [ , good_cluster_idx] = intersect(wf.unitIDs, good_clusters);
            [ ,  MUA_cluster_idx] = intersect(wf.unitIDs, MUA_clusters);
            
            function unpack_wfs(g_or_m_idx, goodness)
                figure;
                if strcmp(goodness, 'good')
                    title('good')
                elseif strcmp(goodness, 'MUA')
                    title('MUA')
                end
                
                for i=1:length(g_or_m_idx) % for each good unit (cluster),...
                    wf_idx = find( wf.unitIDs == g_or_m_idx(i)); % getting index of current 'good' unit

                    cluster_id = wf.unitIDs(wf_idx); % getting index of current 'good' unit
                    best_channel = best_channels{best_channels.Cluster_id == cluster_id, 2}; % getting best channel for current unit

                    cur_wf = squeeze(wf.waveForms(wf_idx, :, best_channel+1,:))';

                    subplot(ceil(length(g_or_m_idx) / 4), 4, i);
                    plot(nanmean(cur_wf, 2))
                    title(['Clus ' int2str(cluster_id) ', chan' int2str(best_channel+1)]);
                    new_struct = struct;
                    new_struct.channel = best_channel + 1;
                    new_struct.unit = cluster_id;
                    new_struct.waveform = cur_wf;
                    new_struct.goodness = goodness;
                    spike_waveform_info = [spike_waveform_info new_struct];
                end
            end
            unpack_wfs(good_cluster_idx, 'good'); 
            unpack_wfs(MUA_cluster_idx, 'MUA');
            
        end
        function waveform_analysis(obj)
            fig = figure('Visible','on','Position',[360,200,450,285]);

            values = uicontrol('Style', 'text', 'String', '',...
                        'Tag', 'p2psym',...
                        'Position', [130, 10, 140, 25]);

            type = uicontrol('Style', 'text', 'String', '',...
                        'Tag', 'type',...
                        'Position', [315, 240, 70, 25]);

            status = uicontrol('Style', 'text', 'String', '',...
                        'Tag', 'status',...
                        'Position', [315, 200, 70, 25]);

            progress = uicontrol('Style', 'text', 'String', '',...
                        'Tag', 'progress',...
                        'Position', [315, 170, 70, 25]);

            good    = uicontrol('Style','pushbutton',...
                         'String','Good',...
                         'Position',[315,140,70,25],...
                         'Callback',@good_Callback);

            bad    = uicontrol('Style','pushbutton',...
                         'String','Bad','Position',[315,100,70,25],...
                         'Callback',@bad_Callback);

            next = uicontrol('Style','pushbutton',...
                         'String','next','Position',[315,60,70,25],...
                         'Callback', @next_Callback);

            back = uicontrol('Style','pushbutton',...
                         'String','back','Position',[315,20,70,25],...
                         'Callback', @back_Callback);

            axes('Units','Pixels','Position',[50,60,200,185]); 

            psth = uicontrol('Style','pushbutton',...
                            'String', 'PSTH',...
                            'Position', [50, 10, 70, 25],...
                            'Callback', @psth_Callback);
            
            function good_Callback(source,~) 
            % for if is good
                [s,key] = get_s;
                s.wf_analysis_goodness = 'Good';
                obj.db(key) = s;
                handle = ancestor(source, 'figure');
                status_handle = findobj(handle, 'Tag', 'status');
                set(status_handle, 'String', 'Included: good');                
            end

            function bad_Callback(source,~) 
            % for if is bad
                [s, key] = get_s;
                s.wf_analysis_goodness = 'Bad';
                obj.db(key) = s;
                handle = ancestor(source, 'figure');
                status_handle = findobj(handle, 'Tag', 'status');
                set(status_handle, 'String', 'Included: bad');
            end

            function back_Callback(source,~) 
            % Display contour plot of the currently selected data.
                handle = ancestor(source, 'figure');
                s = increment(-1, handle);
                
                status_handle = findobj(handle, 'Tag', 'status');
                if isfield(s, 'wf_analysis_goodness')
                    set(status_handle, 'String', s.wf_analysis_goodness);
                else
                    set(status_handle, 'String', 'Included: undecided');
                end 
            end

            function next_Callback(source,~) 
                % for getting next
                handle = ancestor(source, 'figure');

                s = increment(1, handle);
                
                status_handle = findobj(handle, 'Tag', 'status');
                if isfield(s, 'wf_analysis_goodness')
                    set(status_handle, 'String', s.wf_analysis_goodness);
                else
                    set(status_handle, 'String', 'Included: undecided');
                end        
            end        
            
            function s = increment(amt, handle)
                obj.count = obj.count+amt;
                
                
                [s, key] = get_s;

                wf = s.spike_waveforms;
                wf_mean = nanmean(wf,2);
                plot(wf_mean, 'LineWidth', 2); hold on;
                plot(wf_mean+nanstd(wf')', 'b');
                plot(wf_mean-nanstd(wf')', 'b');  
                
                type_handle = findobj(handle, 'Tag', 'type');
                set(type_handle, 'String', ['Sorted as:' s.goodness]);
                
                type_handle = findobj(handle, 'Tag', 'progress');
                set(type_handle, 'String', [num2str(obj.count) '/' num2str(length(obj.wf_keys))]);
                                
                [Vmax,Imax] = max(wf_mean);
                [Vmin,Imin] = min(wf_mean);
                s.p2p = abs(Imax - Imin) / s.amplifier_sampling_rate * 1000 ;
                s.sym = abs(Vmax/Vmin);
                plot(Imax, Vmax, 'm*');
                plot(Imin, Vmin, 'm*');
                hold off;
                obj.db(key) = s;
                
                values_handle = findobj(handle, 'Tag', 'p2psym');
                set(values_handle, 'String', ['p2p: ' num2str(s.p2p) '  sym: ' num2str(s.sym)]);
            end
            
            function [s,key] = get_s
                keys = obj.wf_keys;
                key = keys{obj.count};
                s = obj.db(key);
            end
            
            
            function psth_Callback(~,~)
                % will eventually be the callback for the save function
                db = obj.db; keys = obj.wf_keys;
                key = keys{obj.count};
                s = db(key);
                [on, off, wav_files] = obj.get_stim_on_off(key);
                
                % for each class of stim:
                si = s.stim_identities{1};
                usi = unique(si);
                for i = 1:length(usi)
                    % for each stim itself
                    stim_inds = find(strcmp(si,usi{i}));
                    raster_arr = cell(length(stim_inds),1);
                    for j = 1:length(stim_inds)
                        sp_ts = s.spike_timestamps;
                        ind = stim_inds(j);
                        start = on(ind) - 1 * s.amplifier_sampling_rate;
                        stop = off(ind) + 3 * s.amplifier_sampling_rate;
                        
                        raster_arr{j} = ((intersect(...
                            sp_ts(sp_ts > start),...
                            sp_ts(sp_ts < stop))...
                            - start) / s.amplifier_sampling_rate)';                        
                    end
                    
                    figure;
                    subplot(3,1,1);
                    for k = 1:length(wav_files)
                        if contains(usi{i}, wav_files(k).name)
                            curr_wav = wav_files(k);
                        end
                    end
                    
                    spectrogram(... 
                        [nan(2*s.adc_sampling_rate, 1); curr_wav.data; nan(2*s.adc_sampling_rate, 1)],...
                        256, [],[],s.adc_sampling_rate, 'yaxis')
                    colorbar('delete');
                    
                    title([curr_wav.name '    ' s.context '  why on earth does this histogram work']);
                    
                    subplot(3,1,2);
                    [xpoints, ~] = plotSpikeRaster(raster_arr,...
                            'PlotTYpe','vertline', 'XLimForCell', [0 (stop-start)/s.amplifier_sampling_rate]);
                        
                    histo_axes = subplot(3,1,3);    
                    bin = 0.010; % bin size in s
                    histogram(histo_axes, xpoints, (0:bin:(stop-start)/s.amplifier_sampling_rate)); % convert ms to s
                    xlim([0, (stop-start)/s.amplifier_sampling_rate])
                    % construct psth for that stim
                end
            end 

        end
        
        
        
        % put one point per peak, please
        function cleaned_vec = clean_peaks(~, vec)
            i=1;
            while i < (length(vec))
                if vec(i+1) - vec(i)  < 5
                    vec(i) = [];
                end
                i = i + 1;
            end
            cleaned_vec = vec;
        end
        
        
        %     Logic:
        %     1. Load wav files, text files and audio adc files
        %     2. Downsample them to match the recording SR
        %     3. Cross correlate signals
        %     4. Detect peaks within TTL pulses        
        function [timestamps, filecell, adc_sr] = extract_stim_timestamps(obj, curr_dir)
            workingDirectory = curr_dir;
            cd(workingDirectory)
            %% 1
            disp('adc data')
            % adc data
            stim_data_path = fullfile(workingDirectory, 'adc_data.mat');
            S = load(stim_data_path); % S.board_adc, S.adc_sr
            adc_sr = S.adc_sr;
            
            if contains(workingDirectory, 'mdx')
                audioPath = obj.audioPathMDX;
            else
                audioPath = obj.audioPathMDY;
            end

            disp('diff')
            overunder = [diff(S.board_adc(2,:)) 0]; % add 0 at the end bc this cuts one off
            
            disp('cleaning peaks')
            over = find(overunder > 0.75);
            under = find(overunder < -0.75);
            over = obj.clean_peaks(over);
            under = obj.clean_peaks(under);        
            
            % adc audio
            disp('adc audio')
            adc_audio = S.board_adc(1,:);
            norm_adc_audio = adc_audio - mean(adc_audio(1:S.adc_sr));
            
            % text files for stim id's
            text_path = dir('*markers.txt');
            fid = fopen(fullfile(text_path(1).folder, text_path(1).name),'r');
            filecell = textscan(fid, '%s', 'Delimiter', '\n');
            fclose(fid);
            
            if length(over) ~= length(under) || (length(over) ~= length(filecell{1}))
                error('mismatch in number of stim onset and offset times')
            end
            %% 2
            disp('reading audio file')

            wav_files = dir([audioPath '\*.wav']);
            for i = 1:length(wav_files)
                [cur_wav, fs] = audioread(fullfile(wav_files(i).folder,...
                    wav_files(i).name));
                cur_wav_resampled = resample(cur_wav, S.adc_sr, fs);
                wav_files(i).data = cur_wav_resampled(:,1);             
            end
            
            timestamps = over;
%             timestamps = NaN(length(over),1);
%             for i=1:length(over)
%                 curr_wav_path = split(filecell{1}{i},'\');
%                 curr_wav_name = curr_wav_path{end};
%                 s = wav_files(strcmp({wav_files.name},curr_wav_name));
%                 disp(['Cross correlating Stimulus ' int2str(i)])
%                 [co, lag] = xcorr(norm_adc_audio(over(i):under(i)), s.data);
%                 [~,I] = max(abs(co));
%                 lagDiff = lag(I);
% 
%                 timestamps(i) = over(i) + abs(lagDiff); % Note: timestamp is the middle of the stimulus.
%             end
%             plot(adc_audio)
%             hold on;
% 
%             for u = 1:length(timestamps)
%             plot(timestamps(u), 1 ,'r*')
%             end
%             plot(S.board_adc(2,:));
%             hold off;
            % extract_timestamps function ends here
        end
        
        function filter_raw_data(~)
           target_dir = pwd;
           
           [origFiles, origDataPath] = ... % crucial distinction: files vs file (current)
                uigetfile('*.rhd', 'Select an RHD2000 Data File', 'MultiSelect', 'on');
            cd(origDataPath)
            if iscell(origFiles)
                ts_label = origFiles{1};
            else
                ts_label = origFiles;
            end
            %%
            ts_label = ts_label(1:end-4);
            targetFileName = fullfile(target_dir, 'raw_filtered.dat'); 

            % open raw.dat to write
            fid = fopen(targetFileName, 'w'); % open .dat file for writing
            filearray = [];

            %ordering files
            for i = 1:length(origFiles)
                filearray = [filearray dir(char(origFiles{i}))];
            end

            [~, idx] = sort({filearray.date});
            origFiles = origFiles(idx);
            board_adc = [];
            for i=1:length(origFiles)
                [amplifier_data, frequency_parameters, board_adc_data] = read_Intan_RHD2000_file_MML_DJP(...
                    fullfile(filearray(i).folder, filearray(i).name),0);

                disp([num2str(i) ' of ' num2str(length(origFiles))])

                tic
                % only runs once
                if ~exist('a1', 'var')
                    [b1, a1] = butter(3, 300/frequency_parameters.amplifier_sample_rate*2, 'high');
                end
                % filter
                dataRAW = amplifier_data';
                dataRAW = single(dataRAW);

                datr = filter(b1, a1, dataRAW);
                datr = flipud(datr);
                datr = filter(b1, a1, datr);
                datr = flipud(datr);
                datr=datr';
                fwrite(fid, datr(:),'int16'); % append to .dat file
                %     fwrite(fid1a, amplifier_data(:),'int16'); % append to .dat file
                board_adc = [board_adc board_adc_data];
                toc
            end
            fclose(fid);
            cd(target_dir)
            beep
        end   
        
        
    end
end

%%
% keyset = ['11'; '22'; '33']
% valueset = [struct('field', 'value'), struct('field', 'value'), struct('field', 'value')]
% M = containers.Map('KeyType','char','ValueType', 'struct')
% Add key-value pairs to the map.
% M(keyset(1,:)) = valueset(1)