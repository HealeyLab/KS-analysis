function [keys_r] = get_wf_analysis_keys(obj)
% GET_WF_ANALYSIS_KEYS Returns all the keys whose outputs are for waveform
% analysis
%   Detailed explanation goes here
keys_r = {};
keys = obj.db.keys;
for i = 1:length(keys)
    s = obj.db(keys{i});
    if isfield(s, 'wf_analysis_include') && s.wf_analysis_include
        keys_r{length(keys_r)+1}=keys{i};
    end
end

