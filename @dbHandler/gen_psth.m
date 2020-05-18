function [hs, raster_arr, BL_dur] = gen_psth(obj, argin, show, varargin)
%% Automatically generates psth for each stimulus type
if ischar(argin)
    key = argin;
    s = obj.db(key);
elseif isstruct(argin)
    s = argin;
else
    error('multiple dispatch failure');
end

% get stim on off
[on, off, wav_files] = obj.get_stim_on_off(s);

% for each class of stim:
si = s.stim_identities{1};
uniq = unique(si); % remove redundancies
if isempty(varargin)
    usi = uniq(contains(uniq, 'BOS.wav')); % remove non-wav, non-BOS files
else
    usi = uniq(contains(uniq, varargin{1}));
end
assert(~isempty(usi))
% only for stim type

for i = 1:length(usi)
    %% data
    BL_dur = 2;
    stim_inds = find(strcmp(si,usi{i})); % IN SAMPLES
    raster_arr = cell(length(stim_inds),1);

    for j = 1:length(stim_inds)
        sp_ts = s.spike_timestamps ;
        ind = stim_inds(j);

        % Scooping out spikes from around the buffered stimulus
        start = on(ind) - BL_dur * s.amplifier_sampling_rate;
        stop = off(ind) + BL_dur * s.amplifier_sampling_rate;

        % range in samples, convert back to seconds once returned.
        raster_arr{j} = obj.sliceTS(sp_ts, [start stop]) / s.amplifier_sampling_rate;
    end
    
    %%
    if show
        % formatting output
        for k = 1:length(wav_files)
            if contains(usi{i}, wav_files(k).name)
                curr_wav = wav_files(k);
            end
        end
        figure;
        hs = obj.generatePSTH(raster_arr, s.amplifier_sampling_rate,...
            [zeros(2*s.adc_sampling_rate, 1); curr_wav.data; zeros(2*s.adc_sampling_rate, 1)],...
            s.adc_sampling_rate);
    else
        hs=0;
    end

end

end
