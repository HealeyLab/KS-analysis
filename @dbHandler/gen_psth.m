function gen_psth(obj, key)

s = obj.db(key);

[on, off, wav_files] = obj.get_stim_on_off(key);

% for each class of stim:
si = s.stim_identities{1};
uniq = unique(si); % remove redundancies
usi = uniq(contains(uniq, 'wav')); % remove non-wav files
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

    figure('units','normalized','outerposition',[0 0 .9 .8]);
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
    
    key = replace(key, '_', ' ');
    stim_type = replace(curr_wav.name(1:end-4), '_', ' ');
    title({key;...
        [stim_type ' ' s.context]},...
        'FontSize', 20);

    subplot(3,1,2);
    [xpoints, ~] = plotSpikeRaster(raster_arr,...
            'PlotTYpe','vertline', 'XLimForCell', [0 (stop-start)/s.amplifier_sampling_rate]);

    histo_axes = subplot(3,1,3);    
    bin = 0.010; % bin size in s
    histogram(histo_axes, xpoints, (0:bin:(stop-start)/s.amplifier_sampling_rate)); % convert ms to s
    xlim([0, (stop-start)/s.amplifier_sampling_rate])
%     psth_filename = ['C:\Users\danpo\Documents\dbh_imgs\'...
%         key '_psth_' curr_wav.name(1:end-4) '_' s.context];
%     saveas(fig, psth_filename, 'svg')
% %             saveas(fig, [psth_filename '.pdf'])
%     close(fig)
end

end
