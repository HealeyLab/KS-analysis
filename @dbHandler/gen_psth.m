function [hs, raster_arr] = gen_psth(obj, key)
%% Automatically generates psth for each stimulus type
s = obj.db(key);

[on, off, wav_files] = obj.get_stim_on_off(key);

% for each class of stim:
si = s.stim_identities{1};
uniq = unique(si); % remove redundancies
usi = uniq(contains(uniq, 'BOS.wav')); % remove non-wav files
% only for stim type

for i = 1:length(usi)
    
    %% data
    BL_dur = 2;
    stim_inds = find(strcmp(si,usi{i})); % IN SAMPLES
    raster_arr = cell(length(stim_inds),1);

    for j = 1:length(stim_inds)
        sp_ts = s.spike_timestamps ;
        ind = stim_inds(j);

        start = on(ind) - BL_dur * s.amplifier_sampling_rate;
        stop = off(ind) + BL_dur * s.amplifier_sampling_rate;

        % range in samples, convert back to seconds once returned.
        raster_arr{j} = obj.sliceTS(sp_ts, [start stop]) / s.amplifier_sampling_rate;

    end    
    %% formatting output
    for k = 1:length(wav_files)
        if contains(usi{i}, wav_files(k).name)
            curr_wav = wav_files(k);
        end
    end
    figure;
    hs = obj.generatePSTH(raster_arr, s.amplifier_sampling_rate,...
        [zeros(2*s.adc_sampling_rate, 1); curr_wav.data; zeros(2*s.adc_sampling_rate, 1)],...
        s.adc_sampling_rate);
    %% title and stuff
    axes(hs.sa)
    key = replace(key, '_', ' ');
    stim_type = replace(curr_wav.name(1:end-4), '_', ' ');
    title({key; [stim_type ' ' s.context]}, 'FontSize', 7);
%     psth_filename = ['C:\Users\danpo\Documents\dbh_imgs\'...
%         key '_psth_' curr_wav.name(1:end-4) '_' s.context];
%     saveas(fig, psth_filename, 'svg')
% %             saveas(fig, [psth_filename '.pdf'])
%     close(fig)
end

end
