function T = extractTable(obj)
%%EXTRACT_PORTMANTEAU Connects the BOS playback trials to the BOS listening
%%trial. Hence portmanteau: playback and vocalization are smooshed together
% To get to this step in analysis:
% run dbHandler.add for each recording of cells you need
% run dbHandler.waveform_analysis to mark adbhll included cells' songs
% run set_playback_syllables for each subject (md*_playback_syllables)
% run dbHandler.waveform_connector for each pair of recordings
% run dbHandler.extract_portmanteau to summarize the data!

% Note: I'm treating successive syllable dbhs as the same syllable for
% analysis. The NCM fires the same for them. I could show that NCM doesn't
% have order selectivity by breaking it up by syllable and showing no
% difference.
latency_db = load('C:\Users\danpo\Documents\latency_db.mat');
latency_db = latency_db.latency_db;
field_cell = fields(latency_db);

%% Using a table to aggrgate data
% Columns are Z score, baseline, latency, rows are the pairs of cells, and
% each cell will be a tuple of (playback, song)
% I'll add rows and delete the ones that contain unacceptable values like
% nan or inf

T = table('Size', [1 12], 'VariableTypes',...
    {'string','logical','double','double','double','double','double'     ,'double'      ,'double','double','double','double'},...
    'VariableNames',...
    {'key'   ,'p2p',    'Z_pb' , 'Z_son' ,'FR_pb' ,'FR_son','spkCount_pb','spkCount_son','BL_pb' ,'BL_son','lat_pb','lat_son'});
%     1        2         3        4        5        6        7             8              9        10       11       12                      

%%
% for each field, comprising a pair of keys,...
for i = 1:length(field_cell)
    % convert back to the two fields
    orig_field = field_cell{i};
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
    
    % Match a key in the database
    key1 = char(obj.get_keys(key1));
    subject = obj.get_subject_id(key1);
    key2 = char(obj.get_keys(key2));
    
    % each of these is a struct of syllable objects
    BOS_latencies = latency_db.(orig_field).BOS; 
    SONG_latencies = latency_db.(orig_field).song; 
    
    % the reference for which syllables are which comes from BOS, not song
    all_syllables = fields(BOS_latencies);
    
    % initializing holder variables to be put in table
    Z         = nan(length(all_syllables),1);
    FR        = Z;
    spkcount  = Z;
    latencies = Z;
    
    %% Initialize new row of table
    % { 1     2      3       4        5        6          7          8          9        10       11        12}
    % {'key','p2p', 'Z pb', 'Z son', 'FR pb', 'FR son',  'spkCount','spkCount','BL pb', 'BL son','Lat pb', 'Lat son'}
    newrow = cell(1,12);
    
    %% Playback
    for type_of_syllable = 1:length(all_syllables)
        
        % 'A', 'B', ... 'intro'
        syllable_sTs = BOS_latencies.(all_syllables{type_of_syllable}); % in samples
        
        % Latency
        latency_arr = NaN(length(syllable_sTs),1);
        % Z
        S = zeros(length(syllable_sTs),1);
        numspikes = S; % DJP
        % get array of canonical Syllable objects for BOS of current subject.
        syl_arr = obj.get_playback_syllables(subject).BOS;
        % filter for only the current type of syllable
        syl_arr = syl_arr(contains({syl_arr.id},all_syllables{type_of_syllable}));
        % sort by time
        [~,sort_ind] = sort(arrayfun(@(v)v.window_s(1), syl_arr)); % sort syl_arr according to order or window_s
        syl_arr = syl_arr(sort_ind);
        
        % get S and average latency
        for k = 1:length(syllable_sTs)
            curr_row = syllable_sTs{k};
            % zero truncation, for if some values are below zero due to the cross-correlation. We are going to skip this for now.
%             curr_row = curr_row(curr_row > 0); 
            
            % Z
            if ~isempty(curr_row)
                % S is short for stimulus-evoked activity. It is in FR,
                % which is going to be relativley high.
                % S is the count
                numspikes(k) = length(curr_row);
                curr_syl_ind = mod(k+length(syl_arr)-1, length(syl_arr))+1;
                S(k) = numspikes(k)/diff(syl_arr(curr_syl_ind).window_s); % window is in secs
                % else is already zero
            end
            
            % Latency
            curr_row=curr_row(curr_row>0);
            if ~isempty(curr_row)
                latency_arr(k) = curr_row(1);
                % else is already NaN
            end
        end
        
        % Z scores
        % say there are 300 rows. There are syllables A, B, C, in the order A, A,
        % B, C, B, C. That means you need to find the order of the syllables, and 
        % for that order, plus all the multples of that order, get the zero-truncated
        % window (the spikes have been zero-truncated), and get the FR.
        % FR/Z
        B = obj.get_baseline(key1, 'BOS.wav');
        B = reshape(repmat(B',length(syl_arr),1), [length(syllable_sTs),1]);
        

        % matching sample size
        B = B(1:end);
        S = S(1:end,1);
        covariance = nancov(S(:,1),B);
        Z(type_of_syllable) = (nanmean(S) - nanmean(B)) / sqrt(nanvar(S) + nanvar(B) - 2 * covariance(1,2));
        FR(type_of_syllable) = nanmean(S);
        spkcount(type_of_syllable) = mean(numspikes);
        latencies(type_of_syllable) = nanmean(latency_arr);
    end
%{'key'   ,'p2p',    'Z_pb' , 'Z_son' ,'FR_pb' ,'FR_son','spkCount_pb','spkCount_son','BL_pb' ,'BL_son','lat_pb','lat_son'});

    newrow{3} = Z; % Z-scores
    newrow{5} = FR; % firing rate
    newrow{7} = spkcount; % spike count
    newrow{9} = nanmean(B); % baseline  
    newrow{11} = latencies;
    %% SONG
    % for each type of syllable:
    
    % reinitialize holder variables
    Z         = nan(length(all_syllables),1);
    FR        = Z;
    spkcount  = Z;
    latencies = Z;
    
    for type_of_syllable = 1:length(all_syllables)
        if ~isfield(SONG_latencies, (all_syllables{type_of_syllable}))
            % Sometimes I didn't find any of the current syllable type
            continue
        end
        syllable_sTs_win = SONG_latencies.(all_syllables{type_of_syllable}); %'A', 'B', ... 'intro'
        % Latency
        latency_arr = NaN(length(syllable_sTs_win),1);
        % Z
        S = zeros(length(syllable_sTs_win),1);
        numspikes = S;
        % for each syllable of that type:
        limit = numel(SONG_latencies.(all_syllables{type_of_syllable}).sTs);
        
        for k = 1:limit
            curr_row = syllable_sTs_win.sTs{k};
            % zero truncation. We are going to skip this for now.
%             curr_row = curr_row(curr_row > 0); 
            % Z
            if ~isempty(curr_row) 
                numspikes(k) = length(curr_row);
                S(k) = numspikes(k)/diff(syllable_sTs_win.window_s{k});
                % else is already zero
            end
            
            % Latency
            curr_row = curr_row(curr_row>0);
            if ~isempty(curr_row)
                latency_arr(k) = curr_row(1); % else is already NaN
            end
        end
        
        % Z scores
        % say there are 300 rows. There are syllables A, B, C, in the
        % order A, A, B, C, B, C. That means you need to find the order of
        % the syllables, and for that order, plus all the multples of that
        % order, get the zero-truncated window (the spikes have been
        % zero-truncated), and get the FR. Jesus.
        % FR/Z
        B = obj.db(key2).spike_timestamps;
        B = length(B(B < 2*obj.db(key2).amplifier_sampling_rate)) / 2; % FR for first two seconds
        
        B = repmat(B',length(S),1);
        
        % filter        
        % in case S and B are just zero, increase dimensionality so that 
        % var(1,2) is in bounds
        covariance = nancov(S,B) .* [1 1;1 1]; 
        Z(type_of_syllable) = (nanmean(S) - nanmean(B)) / sqrt(nanvar(S) + nanvar(B) - 2 * covariance(1,2));
        FR(type_of_syllable) = nanmean(S);
        spkcount(type_of_syllable) = mean(numspikes);
        latencies(type_of_syllable) = nanmean(latency_arr);

    end
% {'key'   ,'p2p',    'Z_pb' , 'Z_son' ,'FR_pb' ,'FR_son','spkCount_pb','spkCount_son','BL_pb' ,'BL_son','lat_pb','lat_son'});
    newrow{4}  = Z;
    newrow{6}  = FR;
    newrow{8}  = spkcount;
    newrow{10} = nanmean(B);
    newrow{12} = latencies;
    %% update table
    newrow{1} = orig_field;
    getfirst = @(arr) logical(arr(1));
    newrow{2} = getfirst(obj.get_color(obj.db(key1))); % 1 is narrow, 0 is broad
    T = [T;newrow];
end
%% Clean data
% remove null row
T(1,:)=[]; % removes the first row, which is empty
% make all infs nans
    function cleaned = cleanInf(dirty)
        %% CLEANINF turns infs to nans
        % 
        if iscell(dirty)
            cleaned = cell2mat(dirty);
        end
        cleaned(~isfinite(cleaned)) = nan;
        cleaned = {cleaned};
    end
% clean Z scores
for ii = 1:height(T)
    
    T(ii,3).Z_pb = cleanInf(T(ii,3).Z_pb);
    T(ii,4).Z_son = cleanInf(T(ii,4).Z_son);
%     T(ii,4).BL_pb = 
%     T(ii,5).BL_son = 
%     T(ii,6).lat_pb
%     T(ii,7).lat.son
end
%%
% Z-score
% subplot(131)
vals = [cellfun(@nanmean, T.Z_pb) cellfun(@nanmean, T.Z_son)];

% Note: pointplot eliminates "rungs" in the data that are dirty. That's why
% theres a lot less data in the raw point plots than in the beeswarms.

%
    function [NSpb, NSson, BSpb, BSson] = get_from_table(subplot_ax_ind, T)
        %% GET_FROM_TABLE Breaks out output from the table into four conditions.
        % For subplot_ax_ind = 1, you are getting the third and fourth
        % index, for " = 2, fifth and sixth etc. In other words, for
        % values 1, 2, or 3, you're getting data for Z score, firing rate,
        % or spike number, respectively.
        
        NSpb  = T(T.p2p,subplot_ax_ind*2+1).Variables; % .Variables converts to a cell array
        NSson = T(T.p2p,subplot_ax_ind*2+2).Variables; 
        BSpb  = T(~T.p2p,subplot_ax_ind*2+1).Variables; 
        BSson = T(~T.p2p,subplot_ax_ind*2+2).Variables;
        % fancy way to assign multiple values on a single line
        if iscell(NSpb)
            vc = @(x) vertcat(x{:}); % just flattens the cell array into a vector that can be plotting below
            [NSpb, NSson, BSpb, BSson] = deal(vc(NSpb), vc(NSson), vc(BSpb), vc(BSson));
        end
    end

% Comparison of the population, which is slightly dirtier (not all match up
% nicely between conditions)
% cite: Holger Hoffmann (2020). Violin Plot (https://www.mathworks.com/matlabcentral/fileexchange/45134-violin-plot), MATLAB Central File Exchange. Retrieved March 16, 2020.
figure; hold on
for subplot_axis_ind=1:3
    
    subplot(1,3,subplot_axis_ind)
    [NSpb, NSson, BSpb, BSson] = get_from_table(subplot_axis_ind, T);
       
    plotSpread({NSpb NSson BSpb BSson},...
        'xNames',{'NS' 'BS' 'NS' 'BS'},'distributionMarkers', {'o','o','v','v'},...
        'distributionColors', {obj.NN obj.NN obj.BB obj.BB});

    violin({NSpb NSson BSpb BSson},...
        'xlabel', {'NS', 'NS', 'BS', 'BS'},...
        'facecolor', [obj.NN; obj.NN;obj.BB; obj.BB],'edgecolor', 'none');
    
    % true means NS false means BS
    pointPlot([NSpb NSson ], true(length(NSpb), 1), {'pback' 'song'},...
        { 'NS' 'NS'},{obj.NN obj.NN}, 'xtick', [1,2])
    pointPlot([BSpb BSson ], false(length(BSpb), 1), {'pback' 'song'},...
        {'BS' 'BS'},{obj.BB obj.BB}, 'xtick',[3,4])
    
    
    legend('off')
    xlim([0,5]); % flanking 1 and 4
%     set(gca, 'XTick', xtick)
    set(gca, 'XTickLabel', {'pback' 'song' 'pback' 'song'})
    if subplot_axis_ind == 1
        title('Z-scores')
    elseif subplot_axis_ind == 2
        title('FR')
    elseif subplot_axis_ind == 3
        title('Spike Count')
    end
end
end