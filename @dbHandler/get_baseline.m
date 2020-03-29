function [avg_fr] = get_baseline(obj,key, pattern)
%GET_BASELINE Returns the FR for the baseline for playback trials
%   Detailed explanation goes here
sp_ts = obj.db(key).spike_timestamps;
stim_ids = obj.db(key).stim_identities{1};
[on, off, ~] = obj.get_stim_on_off(key);
indx = find(contains(stim_ids, pattern));
on = on(indx); off = off(indx);
adc_sr = obj.db(key).amplifier_sampling_rate;
window = 3 * adc_sr;
buffer  = 1 * adc_sr;
net_window = window - buffer;

baseline_sTs = cell(length(on),1);
for i = 1:length(on)
    baseline_sTs{i} = length(intersect(...
        sp_ts(sp_ts >= (on(i) - window)),...
        sp_ts(sp_ts < on(i) - buffer))) / (net_window / adc_sr);
end
avg_fr = cellfun(@(c)c(1), baseline_sTs);
end

