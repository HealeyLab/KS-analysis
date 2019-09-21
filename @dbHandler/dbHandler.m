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
% Finally, put everything into a database
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
        audioPaths = {...
            'C:\Users\danpo\Documents\MATLAB\ephysSuite\zf son mdx',...
            'C:\Users\danpo\Documents\MATLAB\ephysSuite\zf son mdy',...
            'C:\Users\danpo\Documents\MATLAB\ephysSuite\zf son mda',...
            'C:\Users\danpo\Documents\MATLAB\ephysSuite\zf son mdb'};
        
        db = containers.Map; % initialized as none, basically
        count = 0;

        wf_keys = {} % a list of keys to use for wf_analysis
        
        BB = [0 63 92] / 255;
        BN = [87 82 126] / 255; 
        NB = [157 99 136] / 255;
        NN = [209 127 127] / 255;
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
        
        function get_song_activity(obj)
            % NOTE: to go fast, toggle the fields with the TAB key.
            % add components

            hs = addcomponents;
            function hs = addcomponents                                             %   [from_left from_bottom width height]
                hs.fig = figure('Visible', 'on', 'Tag', 'fig', 'Units', 'Normalized', 'Position', [0.01 0.2 .95 .7]);
                % open files
                hs.open = uicontrol(hs.fig,...
                    'Style', 'pushbutton','Units', 'Normalized',...
                    'Position',[0.01, 0.01 0.06 0.04],...
                    'String','Open',...
                    'UserData', struct('board_adc',[],'Fs',NaN, 'key',''),...
                    'Callback', @open, 'Tag', 'open');

                % same line
                hs.subjectString = uicontrol(hs.fig,...
                    'Style', 'pushbutton', 'Units', 'Normalized',...
                    'Position',[0.05 0.01 0.3 0.04],... % same line
                    'Style', 'text', 'Units', 'Normalized',...
                    'String','Subject name' ); % ("mda", "mdb", ...)'

                hs.subjectEdit= uicontrol(hs.fig,...
                    'Units', 'Normalized',...
                    'Position',[0.25 0.01 0.05  0.04],... % same line
                    'Style', 'edit',...
                    'String', 'md?', 'Tag', 'subject');

                % Set the save path
                hs.savePathStr= uicontrol(hs.fig,...
                    'Units', 'Normalized', 'Style', 'text',...
                    'Position',[0.33 0.01 0.09 0.04],... 
                    'String','Save path: ');
                hs.savePathEdit = uicontrol(hs.fig,...
                    'Units', 'Normalized', 'Style', 'edit',...
                    'Position',[0.4 0.01 0.3 0.04],...
                    'String', 'C:\Users\danpo\Documents\MATLAB\ephysSuite\',...
                    'Tag', 'savePath');

                % save current selection
                hs.save = uicontrol(hs.fig,...
                    'Units', 'Normalized', 'Style', 'pushbutton',...
                    'Position',[0.9 0.01 0.05  0.04],...
                    'String','Save','Callback', @save,...
                    'Tag', 'save');        

                %         [from_left from_bottom width height]

                hs.SpectroWindow = axes('Units','Normalized', 'Position', [0.04 0.35 0.95 0.75]);
                hs.SpWindow = axes('Units','Normalized', 'Position', [0.04 0.15 0.95 0.15]);
            end

            function save(~, ~)
                function vlines(t)
                    if i == 1
                        col = 'y';
                    elseif i == 2
                        col = 'k';
                    elseif i == 3
                        col = 'b';
                    else
                        col = 'r';
                    end
                    
                    line([floor(t),floor(t)],[0, 1000], 'Color', col, 'LineWidth', 1);
                end
                
                cla(hs.SpWindow);
                
                OpenUD = hs.open.UserData;
                keys = OpenUD.keys;
                for i = 1:length(OpenUD)
                    hold on;
                    key = OpenUD(i).keys;
                    entry = obj.db(key);

                    % get the time limits in samples
                    aFs = entry.amplifier_sampling_rate;
                    axes(hs.SpectroWindow);
                    xl = xlim * 60 * aFs;
                    axes(hs.SpWindow);
                    
                    sTs = entry.spike_timestamps; % query value
                    sTs = sTs(sTs > xl(1) & sTs < xl(2)); % only get within this window
                    sTs = sTs; % - xl(1); % normalize to the beginning of the thing

                    for j = 1:length(sTs)
                        vlines(sTs(j))
                    end
                end
%                 Fs = OpenUD.Fs;
%                 board_adc = OpenUD.board_adc;                
%                 xl = xlim * 60 * Fs; % to audio samples
%                 y = board_adc(end, xl(1):xl(2));
%                 y = y-mean(y);
%                 plot(y)
%                 audiowrite('BOS.wav', y, Fs); 
            end

            function open(hObject, ~)
                % remove all keys not in key family

                keyOrig = input('gimme da key');
                key_fam = strsplit(keyOrig, '&');
                key_fam = key_fam{1};

                keys = obj.db.keys;
                
                fin = false;
                i = 1;
                while ~fin
                    whole_key = keys{i};
                    key_cell = strsplit(whole_key, '&');
                    key_str = key_cell{1};
                    if ~strcmp(key_str, key_fam) 
                        keys = setdiff(keys, {whole_key}); % winnow down
                        i = i - 1;
                    end
                    
                    i = i + 1;
                    if i > length(keys)
                        fin = true;
                    end
                end
                entry = obj.db(keyOrig);
                mic = entry.microphone;
                adc_sr = entry.adc_sampling_rate;
                hObject.UserData = struct('board_adc', mic, 'Fs', adc_sr, 'keys', keys); % frequency_parameters.board_adc_sample_rate
                
                axes(hs.SpectroWindow);
                spectrogram( mic(end,:), 256, [],[], adc_sr, 'yaxis')
                colorbar('delete');
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

        %% Make an entry for each unit in the recording
        % Run per recording
        function add(obj) % to be called as dbHandler.add
            fclose('all');            
            workingDirectory = pwd;
            has_stim_markers = ~isempty(dir('*markers.txt'));
            if has_stim_markers
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
            elseif contains(workingDirectory, 'mda')
                hemisphere = 'L';
            else
                hemisphere = '';
                while ~(strcmp(hemisphere, 'mdy') || strcmp(hemisphere, 'mdx'))
                    hemisphere = input('which hemisphere? (L/R)\n');
                end 
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
            
            % for each element in spike_waveform_info, which has info for
            % all units, 'good' and 'MUA'
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
                
                % for microphone trace
                stim_data_path = fullfile(workingDirectory, 'adc_data.mat');
                S = load(stim_data_path); % S.board_adc, S.adc_sr
                s.microphone = S.board_adc(3,:);
                s.adc_sampling_rate = S.adc_sr;

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
            gwfparams.spikeTimes = readNPY([gwfparams.dataDir '\spike_times.npy']); % ceil(spike_timestamps);
            gwfparams.spikeClusters = readNPY([gwfparams.dataDir '\spike_clusters.npy']); % sp.clu;
            %% Beata's stuff
            gwfparams.chanMap = readNPY([gwfparams.dataDir '\channel_map.npy']); % this is important in esp if you go rid of files. 
            gwfparams.cluster_quality=tdfread([gwfparams.dataDir '\cluster_info.tsv']); % has a bunch of info here. 
            best_channels = (gwfparams.cluster_quality.channel);  
            
            %TODO: add MUA functionality
            good_clusters =gwfparams.cluster_quality.id(...
                find(gwfparams.cluster_quality.group(:,1)=='g')); % | gwfparams.cluster_quality.group(:,1)=='m'));
            gwfparams.good_clusters = good_clusters;
            wf = getWaveForms_BK(gwfparams);
            
            %%
%             wf = getWaveForms(gwfparams);
            %%
            input('when you run the export best channels script, click enter')
            %%
            % CLUSTER GROUP DATA
            if ~isempty(dir('cluster_groups.csv'))
                cluster_quality = readtable('cluster_groups.csv');
            else
                cluster_quality_s = tdfread('cluster_group.tsv', '\t');
                
                Cluster_id = cluster_quality_s.cluster_id;
                % have to put this into a cell array so it's converted to table properly
                group = cell(length(cluster_quality_s.group),1);
                for j=1:length(cluster_quality_s.group)
                    g = cluster_quality_s.group(j,:);
                    % remove spaces
                    g(strfind(g,' ')) = [];
                    group{j} = g;
                end
                
                cluster_quality = table(Cluster_id, group);
            end
            
            % BEST CHANNEL DATA
            %% phy1
            spike_waveform_info = []; % empty struct that will increase in size 
            function unpack_wfs(g_or_m_idx, goodness, spike_waveform_info, wf, best_channels)
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
            if ~isempty(dir('best_channels.csv'))
                best_channels = readtable('best_channels.csv'); % REPLACE THIS

                good_clusters = cluster_quality(strcmp(cluster_quality.group, 'good'),1); % all(cluster_quality.group=='good ', 2); %
                MUA_clusters  = cluster_quality(strcmp(cluster_quality.group, 'mua'),1);
                good_clusters = table2array(good_clusters);
                MUA_clusters  = table2array(MUA_clusters);
                
                [ , good_cluster_idx] = intersect(wf.unitIDs, good_clusters);
                [ ,  MUA_cluster_idx] = intersect(wf.unitIDs, MUA_clusters);

                unpack_wfs(good_cluster_idx, 'good', spike_waveform_info, wf, best_channels); 
                unpack_wfs(MUA_cluster_idx, 'MUA', spike_waveform_info, wf, best_channels);
            %% phy2
            elseif ~isempty(dir('cluster_info.tsv'))
                tsv=tdfread('cluster_info.tsv');
                Cluster_idPhy2 = tsv.id;
                Best_channelPhy2 = tsv.channel;
                goodnessPhy2 = tsv.group;
                
                tsvData = table(Cluster_idPhy2, Best_channelPhy2, goodnessPhy2);
                rows = all(tsvData.goodnessPhy2 == 'good ',2);
                tsvDataGood = tsvData(rows,:);
                rows = all(tsvData.goodnessPhy2 == 'mua  ',2);
                tsvDataMUA = tsvData(rows,:);
                
                for j = 1:height(tsvDataGood)
                    nStrc = struct;
                    nStrc.channel = tsvDataGood(j, :).Best_channelPhy2 - 7; % CUSTOM OFFSET, hardcoded
                    nStrc.unit = tsvDataGood(j, :).Cluster_idPhy2;
                    wfInd = find(wf.unitIDs == nStrc.unit);
                    nStrc.waveform = squeeze(wf.waveForms(wfInd, :, nStrc.channel,:))';
                    nStrc.goodness = tsvDataGood(j, :).goodnessPhy2;
                    spike_waveform_info = [spike_waveform_info nStrc];
                    
                    figure; plot(nanmean(nStrc.waveform, 2))
                end
            else
                error('no bestchannel data')
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
    end
end

%%
% keyset = ['11'; '22'; '33']
% valueset = [struct('field', 'value'), struct('field', 'value'), struct('field', 'value')]
% M = containers.Map('KeyType','char','ValueType', 'struct')
% Add key-value pairs to the map.
% M(keyset(1,:)) = valueset(1)