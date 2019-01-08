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

%%
classdef dbHandler
    properties
        dbPath = 'C:\Users\danpo\Documents\MATLAB\db.mat';
        audioPathMDX = 'C:\Users\danpo\Documents\MATLAB\ephysSuite\zf son mdx';
        audioPathMDY = 'C:\Users\danpo\Documents\MATLAB\ephysSuite\zf son mdy';
        db = containers.Map; % initialized as none, basically
    end
    methods % Static means that methods of this class are not
       % associated with instances of this class.
        function obj = dbHandler()
            S = load(obj.dbPath); 
            obj.db = S.db;
        end

        % Make a key for the Map to point to the unit info
        function key = keyhash(~, workingDirectory, unit, channel)
            label = split(workingDirectory, '\');
            label = label{end};
            % key structure is dir&unit#&channel#, can split on '&'
            key = [label '&' num2str(unit) '&' num2str(channel)];
        end

        % Run per recording
        % Make an entry for each unit in the recording
        function add(obj) % to be called as dbHandler.add(pwd)
            workingDirectory = pwd;
            
            [stim_timestamps, stim_identities, adc_sr] = obj.extract_stim_timestamps();
            [spike_waveform_info, spike_timestamps, sp] = obj.getWaveFormsDriver();
            %%
            % get depth
            depth = input('depth of probe (mm)\n');
            % get context
            context = '';
            while ~(strcmp(context, 'habituation') || strcmp(context, 'random')...
                    || strcmp(context, 't1') || strcmp(context, 't2')...
                    || strcmp(context, 't2') || strcmp(context, 't4'))
                context = input('habituation, random, or t1/t2/t3/t4?\n','s');
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
            
            % for each element in spike_waveform_info, which has info for
            % all units, 'good' and 'MUA'
            for i = 1:length(spike_waveform_info)
                key = obj.keyhash(workingDirectory, spike_waveform_info(i).unit,...
                    spike_waveform_info(i).channel);
                
                obj.db(key) = struct;
                s = obj.db(key);
                s.unit = spike_waveform_info(i).unit;
                s.channel = spike_waveform_info(i).channel;
                s.folder = workingDirectory; % encapsulates day and such
                s.goodness = spike_waveform_info(i).goodness; % 'MUA' or 'good'
                 
                % need each unit's timestamps, not all of them
                s.spike_timestamps = spike_timestamps(sp.clu == spike_waveform_info(i).unit);
                s.amplifier_sampling_rate = sp.sample_rate;
                % actual waveform storage
                s.spike_waveforms = spike_waveform_info(i).goodness;
                
                % complicated ones
                s.depth = depth;
                s.context = context;
                s.hemisphere = hemisphere;
                % adc stuff
                s.stim_timestamps = stim_timestamps;
                s.stim_identities = stim_identities;
                s.adc_sampling_rate = adc_sr;
            end
            save(obj.dbPath, 'db','-v7.3');

        end
        
        function [spike_waveform_info, spike_timestamps, sp] = getWaveFormsDriver(~)
            
            FRAC = 0.10; % fraction of all spikes to sample from
            dataDir = pwd;
            sp = loadKSdir(dataDir);
            
            gwfparams.sr = sp.sample_rate;
            gwfparams.dataDir = dataDir;
            gwfparams.fileName = 'raw.dat';
            gwfparams.dataType = sp.dtype;
            gwfparams.nCh = sp.n_channels_dat;
            gwfparams.wfWin = [-(0.001*gwfparams.sr) 0.002*gwfparams.sr];  % Number of samples before and after spiketime to include in waveform
            gwfparams.nWf = floor(length(sp.st) * FRAC);
            
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
            spike_waveform_info = struct; % empty struct that will increase in size 
            [ , good_cluster_idx] = intersect(wf.unitIDs, good_clusters);
            [ ,  MUA_cluster_idx] = intersect(wf.unitIDs, MUA_clusters);
            
            function unpack_wfs(g_or_m_idx, goodness)
                figure;
                for i=1:length(g_or_m_idx) % for each good unit (cluster),...
                    wf_idx = find( wf.unitIDs == g_or_m_idx(i)); % getting index of current 'good' unit

                    cluster_id = wf.unitIDs(wf_idx); % getting index of current 'good' unit
                    best_channel = best_channels{best_channels.Cluster_id == cluster_id, 2}; % getting best channel for current unit

                    cur_wf = wf.waveFormsMean(wf_idx, best_channel+1,:);
                    cur_wf = squeeze(cur_wf);

                    subplot(ceil(length(g_or_m_idx) / 4), 4, i);
                    plot(cur_wf)
                    title(['Cluster ' int2str(cluster_id) ', channel' int2str(best_channel+1)]);
                    spike_waveform_info(i).channel = best_channel + 1;
                    spike_waveform_info(i).unit = cluster_id;
                    spike_waveform_info(i).waveform = cur_wf;
                    spike_waveform_info(i).goodness = goodness;
                end
            end
            unpack_wfs(good_cluster_idx, 'good');
            unpack_wfs(MUA_cluster_idx, 'MUA');
            
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
        function [timestamps, filecell, adc_sr] = extract_stim_timestamps(obj)
            workingDirectory = pwd;
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
            over = find(overunder > 0.5);
            under = find(overunder < -0.5);
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
            
            timestamps = NaN(length(over),1);
            for i=1:length(over)
                curr_wav_path = split(filecell{1}{i},'\');
                curr_wav_name = curr_wav_path{end};
                s = wav_files(strcmp({wav_files.name},curr_wav_name));
                disp(['Cross correlating Stimulus ' int2str(i)])
                [co, lag] = xcorr(norm_adc_audio(over(i):under(i)), s.data);
                [~,I] = max(abs(co));
                lagDiff = lag(I);
                timestamps(i) = over(i) + abs(lagDiff);
            end
            % extract_timestamps function ends here
        end
        
        function filter_raw_data(~)
           [origFiles, origDataPath] = ... % crucial distinction: files vs file (current)
                uigetfile('*.rhd', 'Select an RHD2000 Data File', 'MultiSelect', 'on');
            cd origDataPath
            if iscell(origFiles)
                ts_label = origFiles{1};
            else
                ts_label = origFiles;
            end
            %%
            ts_label = ts_label(1:end-4);
            dataPath = fullfile(origDataPath,[ts_label '_Kilosort']);
            dataFileName = fullfile(dataPath, 'raw.dat'); 

            % open raw.dat to write
            fid = fopen(dataFileName, 'w'); % open .dat file for writing
            filearray = [];

            %ordering files
            for i = 1:length(origFiles)
                filearray = [filearray dir(char(origFiles{i}))];
            end

            [~, idx] = sort({filearray.date});
            origFiles = origFiles(idx);
            for i=1:length(origFiles)
                read_Intan_RHD2000_file_MML(fullfile(filearray(i).folder,...
                    filearray(i).name),0)

                disp([num2str(i) ' of ' num2str(length(origFiles))])

                tic
                % DRY
                amplifier_data_filtered = bandpass(amplifier_data' ,[350, 4900],...
                frequency_parameters.amplifier_sample_rate,...
                'ImpulseResponse', 'fir'); % fir does not shift peak of trough positions

                amplifier_data_filtered = amplifier_data_filtered';
                fwrite(fid, amplifier_data_filtered(:),'int16'); % append to .dat file
                % DRY
                toc
            end
            fclose(fid); 
        end
    end
end

%%
% keyset = ['11'; '22'; '33']
% valueset = [struct('field', 'value'), struct('field', 'value'), struct('field', 'value')]
% M = containers.Map('KeyType','char','ValueType', 'struct')
% Add key-value pairs to the map.
% M(keyset(1,:)) = valueset(1)