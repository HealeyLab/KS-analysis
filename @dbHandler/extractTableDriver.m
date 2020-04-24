function tableOut = extractTableDriver(obj,key,rasterCell, T)
%EXTRACTTABLEDRIVER Automates organization of a table, specifically for pback trials
%{
    Inputs:
        - key: dbHandler key
        - rasterCell: an M x 1 cell with M = number of pb trials. Will
          contain BL buffer
        - T: table to be updated and spat out as tableOut
%}
    [introKey, motifKey] = obj.getMotifKey(obj.get_subject_id(key));
    
    % sorting because motifs are not always put together in order
    motifs = sort(obj.db(motifKey)); intro = sort(obj.db(introKey));
    vocalization = motifs;
    as = obj.get_audio(obj.get_audioPath(key));
    wavData = as(strcmp({as.name}, 'BOS.wav'));
    motifCell = cell(length(rasterCell), size(vocalization,1));
    BL = cell(length(rasterCell), 1); BLRange = [0 1.5];
    % copy BL shape
    postInhibitoryRebound = BL;
    for r = 1:length(rasterCell)
        % Each row is a different playback. Each column is a motif within
        % that playback
        % [number of distinct motifs, t] = size(motifs)
        for m = 1:size(vocalization,1)
            
            % accounting for the 2 seconds added on manually, hardcoded in. 
            % because rasterCell includes 2 seconds buffer
            % I'd like to fix this in the future;
            range = vocalization(m,:) + 2; 
            
            % get spike timestamps for motif
            motifCell{r,m} = obj.sliceTS(rasterCell{r}, range); 
            motifCell{r,m} = motifCell{r,m}'; % needs to be transposed for plotSpikeRaster
        end
        % get BL activity for playback
        BL{r} = obj.sliceTS(rasterCell{r}, BLRange);
        % last motif, move past buffer, and add the motif vector to the
        % [0 1] range
        postInhibitoryRebound{r} = obj.sliceTS(rasterCell{r}, vocalization(end,:) + 2 + [0 1]);
    end
    tableOut = obj.extractTableHelper(T, key, motifCell, vocalization, BL, diff(BLRange), cellfun(@length, postInhibitoryRebound));
    %% add to PDF
    % Above I organized the data to fit into the table in just the right
    % way. Now, I want to vizualize a little bit more, so I'm augmenting
    % the motif and motifCell
    
    % Adding baseline data on either side, but for analysis we are using
    % only the before BL window for analysis
    
    % subtracting one from the first one to signal to genMotifPSTH.m.
    % adding the length of the BLRange to motifs so it's offset properly.
    vocalization = [BLRange + [1 0]; vocalization + diff(BLRange)]; % number of motifs will drive number of PSTHs
    motifCell = [cellfun(@(x) x', BL, 'UniformOutput', false) motifCell];
    %{
    if isempty(dir(['C:\Users\danpo\Documents\dbh_imgs\*' key '.pdf']))
        obj.saveMotifFigs(key, ...
            [zeros(diff(BLRange)*obj.db(key).adc_sampling_rate, 1); wavData.wav; zeros(diff(BLRange)*obj.db(key).adc_sampling_rate, 1)],...
            wavData.fs, obj.db(key).amplifier_sampling_rate, motifCell, motifs);
    end
    %}

end

