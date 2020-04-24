function [T] = extractTableHelper(obj, T_in, key, motifCell, motifs, BL, BL_dur, postInhibitoryRebound)
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
    entry = obj.db(key);    
    
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
    Z = [];
    for i = 1:size(FR,1)
        % To evenly distribute the variance between the song and playback
        % conditions, I going to take a random trial index from the stim FR and
        % the BL FR
        randInd = floor(size(FR, 2) * rand)+1;
        covariance = nancov(FR(i, randInd), BL_FR(randInd)); % ask Luke about this
        if length(covariance) > 1
            covariance = covariance(2,1);
        elseif length(covariance) == 1
            covariance = covariance(1,1);
        else
            assert('covariance is a weird shape')
        end
            
        if ~isempty(motifCell{1})% basically, if there aren't no spikes
            new_Z = (nanmean(FR)-nanmean(BL_FR)) / sqrt(nanvar(FR)+nanvar(BL_FR)-2*covariance);
        else
            new_Z = 0; 
        end
        Z = [Z; new_Z];
    end
    
    win = diff(motifs');    
    newrow = {key, obj.get_p2p(entry), entry.hemisphere, obj.get_sym(entry),...
        Z, FR, spk',0, win', BL_FR, BL_spk, postInhibitoryRebound}; % Z, mean(FR,1)', mean(spk,2),0, win', BL_FR, BL_spk, postInhibitoryRebound};
    
    T = [T_in;newrow];    
end


