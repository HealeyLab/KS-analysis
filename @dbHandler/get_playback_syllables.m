function syls = get_playback_syllables(obj,subject_pattern)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

keys = obj.get_keys('playback_syllables');
inds = cellfun(@(s)contains(s,subject_pattern),keys);
keys = keys(find(inds));
syls = obj.db(char(keys));
end

