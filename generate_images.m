%% generate all images
% dbh = dbHandler();
wf_keys = dbh.wf_keys;
db = dbh.db;

BB = dbh.BB;
BN = dbh.BN;
NB = dbh.NB;
NN = dbh.NN;
%% swatch
% figure
% subplot(2,2,1)
% set(gca, 'Color',BB)
% set(gca, 'YTick', [], 'XTick', [])
% title('BB')
% set(gca, 'fontSize', 30)
% 
% subplot(2,2,2)
% set(gca, 'Color',NN)
% set(gca, 'YTick', [], 'XTick', [])
% title('NN')
% set(gca, 'fontSize', 30)
% 
% subplot(2,2,3)
% set(gca, 'Color',BN)
% set(gca, 'YTick', [], 'XTick', [])
% title('BN')
% set(gca, 'fontSize', 30)
% 
% subplot(2,2,4)
% set(gca, 'Color',NB)
% set(gca, 'YTick', [], 'XTick', [])
% title('NB')
% set(gca, 'fontSize', 30)
%% 
habit_arr = cell(length(dbh.wf_keys),2);
latency_arr = cell(length(dbh.wf_keys),2);
ns = 0;
bs = 0;
for i = 1:length(dbh.wf_keys)
    key = wf_keys{i};

    s = db(key);
    %% waveform
%     wf_filename = ['C:\Users\danpo\Documents\dbh_imgs\' wf_keys{i} '_wf'];
%     fig = figure;
%     wf = s.spike_waveforms;
%     wf_mean = nanmean(wf,2);
%     hold on
%     if s.p2p >= .43
%         color = BB;
%     else
%         color = NN;
%     end
%     plot(wf_mean, 'LineWidth', 1.25, 'Color', color);
%     plot(wf_mean+nanstd(wf')','LineWidth', 1.25, 'Color', color);
%     plot(wf_mean-nanstd(wf')','LineWidth', 1.25, 'Color', color);  
%     set(gca, 'XTick', [])
%     title(['p2p: ' num2str(s.p2p) '            sym: ' num2str(s.sym)])
%     saveas(fig, wf_filename, 'svg')
%     close(fig)
    %% habituation rate (evoked) (now in FR, minus baseline)
    if isfield(s, 'stim_timestamps')
        if s.p2p >= 0.43
            bs = bs + 1;
        else
            ns = ns + 1;
        end 
        [on, off, wav_files] = dbh.get_stim_on_off(key);

        si = s.stim_identities{1};
        usi = unique(si);
        usi_leg = unique(si);
        for j = 1:length(usi_leg)
            edit = usi{j};
            edit = strsplit(edit, '\');
            edit = edit{end};
            edit = edit(1:end-4);
            usi_leg{j} = edit;
        end
        
        % for each class of stim:
        lat_arr = nan(length(usi),4); % latency to first spike for each stim onset
        fr_arr = nan(length(usi), 4); % will overlay response to each of four stims
        for j = 1:length(usi)
            stim_inds = find(strcmp(si,usi{j}));
            % for each stim itself
            for k = 1:length(stim_inds)
                sp_ts = s.spike_timestamps;
                ind = stim_inds(k);
                evoked = length(intersect(...
                    sp_ts(sp_ts >= on(ind)),...
                    sp_ts(sp_ts < off(ind)))) / ((off(ind)-on(ind)) / s.amplifier_sampling_rate);                      
                
                baseline_length_samples = 3 * s.amplifier_sampling_rate;
                baseline_length_buffer  = 1 * s.amplifier_sampling_rate;
                
                baseline = length(intersect(...
                    sp_ts(sp_ts >= (on(ind) - baseline_length_samples)),...
                    sp_ts(sp_ts < on(ind)-baseline_length_buffer))) / (3-1);
                
                fr_arr(j,k) = evoked - baseline; % now in FR!
                % get latency now
                all_greater = sp_ts(sp_ts >= on(ind));
                next_greater = all_greater(1);
                lat_arr(j,k) = (next_greater-on(ind))/s.amplifier_sampling_rate;
            end
        end
       % get legend ready
        habit_arr{i, 1} = fr_arr;
        habit_arr{i, 2} = s.p2p; 
        latency_arr{i,1} = lat_arr;
        latency_arr{i,2} = s.p2p;
        
%         fig = figure('units','normalized','outerposition',[0 0 1 1]);
%         plot(fr_arr', 'LineWidth', 1.25)
%         legend(usi_leg)
%         title('Stimulus-specific adaptation', 'FontSize', 24)
%         xlabel('Stimulus presentation' , 'FontSize', 24)
%         ylabel('Evoked firing rate', 'FontSize', 24)
%         %inset for waveform
%         axes('Position', [.7, .2, .15,.15])
%         wf_mean = nanmean(s.spike_waveforms, 2);
%         wf_std = nanstd(s.spike_waveforms');
%         plot(wf_mean, 'Color', 'k'); hold on;
%         plot(wf_mean+wf_std', 'Color', 'k');
%         plot(wf_mean-wf_std', 'Color', 'k');  
%     
%         title(['p2p: ' num2str(s.p2p)], 'FontSize', 24)
%         %save
%         ssa_filename = ['C:\Users\danpo\Documents\dbh_imgs\'...
%                 key '_ssasub_' s.context];
%         saveas(fig, ssa_filename, 'svg')
%         close(fig)
    end
   
    %% psth
    if 0 %isfield(s, 'stim_timestamps')
        [on, off, wav_files] = dbh.get_stim_on_off(key);

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
                start = on(ind) - 2 * s.amplifier_sampling_rate;
                stop = off(ind) + 2 * s.amplifier_sampling_rate;

                raster_arr{j} = ((intersect(...
                    sp_ts(sp_ts > start),...
                    sp_ts(sp_ts < stop))...
                    - start) / s.amplifier_sampling_rate)';                        
            end

            fig = figure('units','normalized','outerposition',[0 0 1 1]);
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

            title([curr_wav.name(1:end-4) '    ' s.context],'FontSize', 24);

            subplot(3,1,2);
            [xpoints, ~] = plotSpikeRaster(raster_arr,...
                    'PlotTYpe','vertline', 'XLimForCell', [0 (stop-start)/s.amplifier_sampling_rate]);

            histo_axes = subplot(3,1,3);    
            bin = 0.010; % bin size in s
            histogram(histo_axes, xpoints, (0:bin:(stop-start)/s.amplifier_sampling_rate)); % convert ms to s
            xlim([0, (stop-start)/s.amplifier_sampling_rate])
            psth_filename = ['C:\Users\danpo\Documents\dbh_imgs\'...
                key '_psth_' curr_wav.name(1:end-4) '_' s.context];
            saveas(fig, psth_filename, 'svg')
%             saveas(fig, [psth_filename '.pdf'])
            close(fig)
        end
    end
end
 %% analyze habit_arr
habit_arr = habit_arr(~cellfun('isempty', habit_arr));
habit_arr = reshape(habit_arr, [length(habit_arr)/2,2]);

% First, beeswarm, taking average of whole set
ns = [];
bs = [];
for j = 1:length(habit_arr)
    if habit_arr{j,2} >= .43
        bs = [bs mean(mean(habit_arr{j,1}))];
    else
        ns = [ns mean(mean(habit_arr{j,1}))];
    end
end
f = figure;
plotSpread({bs, ns}, 'xNames', {'Broad','Narrow'},...
    'distributionColors', {BB, NN})
set(findall(f, 'type','line'),'markerSize',20)
ylabel('Evoked Firing Rate', 'FontSize', 20)
title('Overall evoked firing rate', 'FontSize', 20)
[h,p] = ttest2(bs, ns);
sigstr = '';
if p<0.001
    sigstr = '***';
elseif p<0.01
    sigstr = '**';
elseif p<0.05
    sigstr = '*';
end
text(1.25, 27, sigstr, 'FontSize', 28)

% next, beeswarm for each stimulus
ns = [];
bs = [];
for j = 1:length(habit_arr)
    ha = habit_arr{j,1};
    stims = nan(4,1);
    for k = 1:4
        if size(ha,1) >= k
            stims(k) = mean(mean(ha(k,:)));
%         else
%             stims(k) = nan;
        end
    end
    
    % bos bosrev con wn
    if habit_arr{j,2} >= .43
        bs = [bs; stims(1) stims(2) stims(3) stims(4)];
    else
        ns = [ns; stims(1) stims(2) stims(3) stims(4)];
    end
end
f = figure;
plotSpread({bs(:,1), ns(:,1), bs(:,2), ns(:,2), bs(:,3), ns(:,3), bs(:,4), ns(:,4)},...
    'xNames', {'BS BOS','NS BOS','BS BOS REV','NS BOS REV','BS CON','NS CON', 'BS WN','NS WN'},...
    'distributionColors', {BB, NN, BB, NN, BB, NN, BB, NN})
set(findall(f, 'type','line'),'markerSize',12);
ylabel('Evoked Firing Rate')
title('Evoked firing rate for each stimulus')
% t-test show
for j = 1:4
    [~,p,~,stats]=ttest2(ns(:,j), bs(:,j))
    sigstr = '';
    if p<0.001
        sigstr = '***';
    elseif p<0.01
        sigstr = '**';
    elseif p<0.05
        sigstr = '*';
    end
    text(2*j-0.75,35,sigstr, 'FontSize', 28)
%         text(2*j-0.5,45,['p=' num2str(p)])
end

% finally, waterfall
bs = nan(4, 80, 88);
ns = bs; % smallen these last dimensions before generating graph
for j = 1:length(habit_arr)
    ha = habit_arr{j,1};
    % stims = nan(numstims, num_experiments, num_presentations_per_experiment)
    stims = nan(4,80);
    for k = 1:4
        if size(ha,1) >= k
            stims(k,:) = [ha(k,:) nan(1,80-size(ha,2))];
        end
    end
    
    % bos bosrev con wn
    if habit_arr{j,2} >= .43
        bs(:,:,j) = stims;
    else
        ns(:,:,j) = stims;
    end
    
end

figure;
title('Stimulus-specific adaptation curves')
usi_leg = {'BOS', 'BOS REV', 'CON', 'WN'};
for i = 1:4
    ns_curve = squeeze(ns(i,:,:));
    bs_curve = squeeze(bs(i,:,:));
    yl = max(max([ns_curve; bs_curve]));
    
    subplot(2,4,i)
    plot(bs_curve, 'Color', BB, 'LineWidth', 1); hold on
    plot(nanmean(bs_curve,2), 'Color', 'r', 'LineWidth', 1.2)
%     figure;waterfall(squeeze(ns(i,:,:))')
    title(usi_leg{i});
    if i == 1
    ylabel('Evoked Firing Rate')
    end
    xlim([1 50])
    ylim([1 yl])
        
    subplot(2,4,i+4)
    plot(ns_curve, 'Color', NN, 'LineWidth', 1); hold on
	plot(nanmean(ns_curve,2), 'Color', 'k', 'LineWidth', 1.2)
%     figure;waterfall(squeeze(bs(i,:,:))')
    title(usi_leg{i});
    xlabel('Stimulus presentation')
    
    if i == 1
        ylabel('Evoked Firing Rate')
    end
    
    xlim([1 50])
    ylim([1 yl])
end
%% Latency analysis
latency_arr = latency_arr(~cellfun('isempty', latency_arr));
latency_arr = reshape(latency_arr, [length(latency_arr)/2,2]);

% First, beeswarm, taking average of whole set
ns = [];
bs = [];
for j = 1:length(latency_arr)
    if latency_arr{j,2} >= .43
        bs = [bs mean(mean(latency_arr{j,1}))];
    else
        ns = [ns mean(mean(latency_arr{j,1}))];
    end
end
f = figure;
% whisker
sorted = [bs ns];
grp = [ones(length(bs),1); 2 * ones(length(ns),1)];
boxplot(sorted, grp)

% beeswarm
plotSpread({bs, ns}, 'xNames', {'Broad','Narrow'},...
    'distributionColors', {BB, NN})
set(findall(f, 'type','line'),'markerSize',18)
ylabel('First spike latency (s)', 'FontSize', 18)
ylim([0 1.2])
title('Overall first spike latency', 'FontSize', 18)

% next, beeswarm for each stimulus
ns = [];
bs = [];
for j = 1:length(latency_arr)
    ha = latency_arr{j,1};
    stims = nan(4,1);
    for k = 1:4
        if size(ha,1) >= k
            stims(k) = mean(mean(ha(k,:)));
%         else
%             stims(k) = nan;
        end
    end
    
    % bos bosrev con wn
    if latency_arr{j,2} >= .43
        bs = [bs; stims(1) stims(2) stims(3) stims(4)];
    else
        ns = [ns; stims(1) stims(2) stims(3) stims(4)];
    end
end
figure
% do whisker
sorted = [bs(:,1); ns(:,1); bs(:,2); ns(:,2); bs(:,3); ns(:,3); bs(:,4); ns(:,4)];
grp = [  1* ones(length(bs(:,1)),1);...
   2* ones(length(ns(:,1)),1);...
   3* ones(length(bs(:,2)),1);...
   4* ones(length(ns(:,2)),1);...
   5* ones(length(bs(:,3)),1);...
   6* ones(length(ns(:,3)),1);...
   7* ones(length(bs(:,4)),1);...
   8* ones(length(ns(:,4)),1)];
boxplot(sorted, grp)

plotSpread({bs(:,1), ns(:,1), bs(:,2), ns(:,2), bs(:,3), ns(:,3), bs(:,4), ns(:,4)},...
    'xNames', {'BS BOS','NS BOS','BS BOS REV','NS BOS REV','BS CON','NS CON', 'BS WN','NS WN'},...
    'distributionColors', {BB, NN, BB, NN, BB, NN, BB, NN})
f = gca;
set(findall(f, 'type','line'),'markerSize',15);
ylabel('First spike latency (s)', 'FontSize', 18)
ylim([0 1.2])
title('First spike latency for each stimulus', 'FontSize', 18)
% independent t-test:
% [h,p]=ttest2(bs(:,1), ns(:,1))
%%
% [BBavg, NNavg, BNavg, NBavg] = dbh.cross_correlograms;
%
%% beeswarms
m = '.';
f=figure;
plotSpread({BBavg(:,1), NNavg(:,1), BNavg(:,1), NBavg(:,1),...
            BBavg(:,2), NNavg(:,2), BNavg(:,2), NBavg(:,2),...
            BBavg(:,3), NNavg(:,3), BNavg(:,3), NBavg(:,3),...
            BBavg(:,4), NNavg(:,4), BNavg(:,4), NBavg(:,4)},....
'xNames', { 'BB','NN','BN','NB',...
            'BB','NN','BN','NB',...
            'BB','NN','BN','NB',...
            'BB','NN','BN','NB'},...
'distributionColors', {BB, NN, BN, NB,  BB, NN, BN,NB,...
                       BB, NN, BN, NB,  BB, NN, BN,NB},...
'distributionMarkers', {m,m,m,m,  m,m,m,m,  m,m,m,m,  m,m,m,m})

% plotSpread({BBavg, NNavg, BNavg, NBavg},....
% 'xNames', {'B/B','N/N','B/N','N/B'}, 'distributionColors', {BB, NN, BN,NB},...
% 'distributionMarkers', {m,m,m,m})

ylabel('Expected value (ms)', 'FontSize', 15)
title('Expected latency across cross-correlations before and after zero', 'FontSize', 15)
set(findall(f,  'type','line'),'markerSize',12)

for i = 4:4:16
    line([i + .5 i + .5], [-50, 50],'Color','k')
end

text(1,25, 'Minimum Left')
text(5,25, 'Maximum Left')
text(9,-25, 'Minimum Right')
text(13,-25, 'Maximum Right')