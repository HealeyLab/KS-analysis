function [introKey, motifKey] = getMotifKey(obj, key)
%GETMOTIFKEY Returns the dbh.db keys for intro and motif windows.
    sn = split(key, 'Kilo');
    introKey = [sn{1} '_intro'];
    motifKey = [sn{1} '_motifs'];
end

