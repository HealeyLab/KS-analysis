function [T] = extractTableHelper(obj, T_in, argin, motifCell, motifs, BL, BL_dur, postInhibitoryRebound)
%EXTRACTTABLEHELPER Helper function for extractTables
%   Does the dirty work of generating a new table
%{
    Input arguments:
        T_in - empty table to add rows to
        key - dbHandler key indicating a unit
        motifCell - m x 2 cell array containing indicies (in seconds) of
            motifs
        motifs - m x 1 cell array containing zero-subtracted spike
            timestamps
        BL - 
        BL_dur - duration of BL
        postInhibitoryRebound - the spikes that occur right after the last
        motif
%}
    if ischar(argin)
        key = argin;
        entry = obj.db(key);    
    elseif isstruct(argin)
        entry = argin;
    else
        error('not a key or an entry')
    end
    
    % convert timestamps to numspks
    spk    = cellfun(@length, motifCell); 
    BL_spk = cellfun(@length, BL);
    
    % Divide numspks by time durations
    FR = (spk ./ diff(motifs'))';
    % if it's a single-line vector, you won't get pair-wise comparison. SO
    % we match the length of BL_spk to that of FR
    BL_spk = ones(size(FR,2),1) .* BL_spk;
    BL_FR = (BL_spk / BL_dur);
    if length(BL_FR) == 1
        BL_FR = ones(length(FR),1) * BL_FR;
    end
    
    % I'm averaging the FR from each motif because I want to pair each
    % playback the baseline to construct a covariance matrix

    win = diff(motifs');
    
    % if there's only one motif
    if length(win) ==  1
        FR = {FR};
        spk = {spk};
        win = {win};
        BL_FR = {BL_FR};
        BL_spk = {BL_spk};
        postInhibitoryRebound = {postInhibitoryRebound};
    end
        
    newrow = {key, obj.get_p2p(entry), entry.hemisphere, obj.get_sym(entry),...
        FR, spk',0, win', BL_FR, BL_spk, postInhibitoryRebound}; % Z, mean(FR,1)', mean(spk,2),0, win', BL_FR, BL_spk, postInhibitoryRebound};
    
    T = [T_in;newrow];    
end


