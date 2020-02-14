function [on, off, wav_files] = get_stim_on_off(obj,key)
    s = obj.db(key);
    audioPath = obj.get_audioPath(s.folder);

    wav_files = dir([audioPath '\*.wav']);
    for i = 1:length(wav_files)
        [cur_wav, fs] = audioread(fullfile(wav_files(i).folder,...
            wav_files(i).name));
        cur_wav_resampled = resample(cur_wav, s.amplifier_sampling_rate, fs);
        wav_files(i).data = cur_wav_resampled(:,1);             
    end
    % aligning sr of stim timestamps
    s.stim_timestamps = s.stim_timestamps / s.adc_sampling_rate * s.amplifier_sampling_rate;
    % wait, stim_timestamps are in a different sampling rate right?
    on = nan(length(s.stim_timestamps),1); off = nan(length(s.stim_timestamps),1);
    for i = 1:length(s.stim_timestamps)
        curr_wav_path = split(s.stim_identities{1}{i},'\');
        curr_wav_name = curr_wav_path{end};
        index = strcmp({wav_files.name},curr_wav_name); % for the non-tone stimuli
        if any(index)
            len = length(wav_files(index).data);
        else
            len = length(tones(200, 1, 30e3));
        end
        on(i) = s.stim_timestamps(i); 
        off(i)= s.stim_timestamps(i) + len;
    end
end
