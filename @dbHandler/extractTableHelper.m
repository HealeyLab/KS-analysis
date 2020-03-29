function [T] = extractTableHelper(obj, T_in, key, motifCell, motifs, sonogram, fs, BL, BLRange)
%EXTRACTTABLEHELPER Helper function for extractTables
%   Does the dirty work of generating a new table
%
    entry = obj.db(key);
    
    
    % hopefully turn the code below into a function.
    figure;
    for m = 1:size(motifs,1)
        % viz
        mot = motifs(m,:) * fs; % converting to samples
        obj.generatePSTH(motifCell(:,m), entry.amplifier_sampling_rate,...
               	 sonogram(mot(1):mot(2)), fs, size(motifs,1), m);
    end
    % aggregating
%     FR(r) = 0
%     covariance = nancov(motifArr{:}, BL);
%     Z(r) = (nanmean(S) - nanmean(BL)) / sqrt(nanvar(S) + nanvar(BL) - 2 * covariance(1,2));
%     spk(r) = 0
    spk = cellfun(@(c) length(c), motifCell); % convert to numspks
    FR = mean(spk ./ diff(motifs'),2); % hey cool this works as column-wise division!
    BL_spk = cellfun(@length, BL);
    BL_FR = BL_spk / diff(BLRange);
    % I'm averaging the FR from each motif because I want to pair each
    % playback the baseline to construct a covariance matrix
    covariance = nancov(FR, BL_FR);
    Z = (nanmean(FR)-nanmean(BL_FR)) / sqrt(nanvar(FR)+nanvar(BL_FR)-covariance(1,2));
    newrow = cell(1,width(T_in));
    newrow = {key, obj.get_p2p(entry), entry.hemisphere, obj.get_sym(entry),...
        Z, FR, spk,0, mean(diff(motifs')), BL_FR, BL_spk};
    
    T = [T_in;newrow];    
end
function drawLines(arr, c)
        hold on;

    for stimInd = 1:size(arr,1)
        line([arr(stimInd,1), arr(stimInd,1)], [0,15e3], 'Color', c);
        line([arr(stimInd,2), arr(stimInd,2)], [0,15e3], 'Color', c);
    end
end

