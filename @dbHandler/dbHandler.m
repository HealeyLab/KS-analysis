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
        dbPath = 'C:\Users\danpo\Documents\db.mat';%'E:\DJP thesis sorting\db.mat';
        audioPaths = {...
            'C:\Users\danpo\Documents\MATLAB\ephysSuite\zf son mdx',...
            'C:\Users\danpo\Documents\MATLAB\ephysSuite\zf son mdy',...
            'C:\Users\danpo\Documents\MATLAB\ephysSuite\zf son mda',...
            'C:\Users\danpo\Documents\MATLAB\ephysSuite\zf son mdb',...
            'C:\Users\danpo\Documents\MATLAB\ephysSuite\zf son mdc',...
            'C:\Users\danpo\Documents\MATLAB\ephysSuite\zf son mdd',...
            'C:\Users\danpo\Documents\MATLAB\ephysSuite\zf son mde'};
        
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
        %%
        function color = get_color(obj, s)
            if obj.get_p2p(s) >= .43
                color = obj.BB;
            else
                color = obj.NN;
            end
        end
        %%
        function figs = generate_figures(obj, key_pattern)
            % return a list of figs that can ust be deleted real easily
            figs = [];
            keys = obj.show_keys(key_pattern);
            for i = 1:length(keys)
                key = keys{i};
                s = obj.db(key);
                
                waveforms = obj.gen_waveforms(key, 1); % 1 means only one
                figs = [figs waveforms];
                if  strcmp(s.context, 'song')
                    son = obj.get_song_activity(key, 1); % 1 means only one
                    figs = [figs son];
                elseif strcmp(s.context, 'habituation')
                    psth = obj.gen_psth(key);
                    strf = obj.gen_strf(key);
                    
                    figs = [figs psth strf];
                end
            end
        end
        %%
        function p2p = get_p2p(~, s)
            wf = s.spike_waveforms;
            wf_mean = nanmean(wf,2);
            [~, min_i] = min(wf_mean);
            [~, max_i] = max(wf_mean);
            p2p = (max_i-min_i) / s.amplifier_sampling_rate * 1000;
        end
        
        function sym = get_sym(~, s)
            wf = s.spike_waveforms;
            wf_mean = nanmean(wf,2);
            [~, min_i] = min(wf_mean);
            [~, max_i] = max(wf_mean);
            sym = max_i / min_i;
        end
        %%
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
                    || strcmp(context, 'song') || strcmp(context, 'other'))
                context = input('habituation, random, song, or other?\n','s');
            end
            
            % get hemisphere
            if contains(workingDirectory, 'mdx')...
                    || contains(workingDirectory, 'mde')
                hemisphere = 'R';
            elseif contains(workingDirectory, 'mdy') ...
                    || contains(workingDirectory, 'mda')...
                    || contains(workingDirectory, 'mdb')...
                    || contains(workingDirectory, 'mdc')...
                    || contains(workingDirectory, 'mdd')
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
            wf = obj.getWaveForms_BK(gwfparams);
            
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
        
        function show_each_entry_size(obj)
            for i = 1:length(obj.db)
                keys = dbh.db.keys;
                s = dbh.db(keys{i});
                disp(keys{i})
                whos('s')
            end
        end
    end
end

%%
% keyset = ['11'; '22'; '33']
% valueset = [struct('field', 'value'), struct('field', 'value'), struct('field', 'value')]
% M = containers.Map('KeyType','char','ValueType', 'struct')
% Add key-value pairs to the map.
% M(keyset(1,:)) = valueset(1)