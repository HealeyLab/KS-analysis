function tableOut = extractTableDriver(obj,key,rasterCell, T, BL_buffer, s, vocalizationType)
%                                          1   2           3  4       5
%EXTRACTTABLEDRIVER Automates organization of a table, specifically for pback trials
%{
    Inputs:
        - key: dbHandler key
        - rasterCell: an M x 1 cell with M = number of pb trials. Will
          contain BL buffer
        - T: table to be updated and spat out as tableOut
        - varargin: the current entry, a struct
%}

    [introKey, motifKey] = obj.getMotifKey(obj.get_subject_id(key));
    % sorting because motifs are not always put together in order
    %{
        Error using containers.Map/subsref
        The specified key is not present in this container
        means: run obj.basicViewer(key)
    %}
    
    % dynamically assign vocalization. I know it's ugly.
    if strcmp(vocalizationType, 'motifs')
        vocalization = sort(obj.db(motifKey));
    elseif strcmp(vocalizationType, 'intro')
        vocalization = sort(obj.db(introKey));
    else
        error('vocalization parameter must be ''motifs'' or ''intro''')
    end
    
    vocalization = vocalization + BL_buffer;
    
    as = obj.get_audio(obj.get_audioPath(key));
    wavData = as(strcmp({as.name}, 'BOS.wav'));
    motifCell = cell(length(rasterCell), size(vocalization,1));
    
    BL_dur = 1; %s
    BL = cell(length(rasterCell), 1); BLRange = [0 BL_dur];
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
            range = vocalization(m,:); 
            
            % get spike timestamps for motif
            motifCell{r,m} = obj.sliceTS(rasterCell{r}, range); 
            motifCell{r,m} = motifCell{r,m}'; % needs to be transposed for plotSpikeRaster
        end
        
        % get BL activity for playback
        BL{r} = obj.sliceTS(rasterCell{r}, BLRange);
        
        % last motif, move past buffer, 
        PIR_dur = 0.25;
        postInhibitoryRebound{r} = obj.sliceTS(rasterCell{r}, vocalization(end,2) + [0 PIR_dur]);
    end
    tableOut = obj.extractTableHelper(T, key, motifCell, vocalization,...
        BL, diff(BLRange), cellfun(@(x) length(x) / PIR_dur, postInhibitoryRebound));

    %{
    %% pdf
    % Above I organized the data to fit into the table in just the right
    % way. Now, I want to vizualize a little bit more, so I'm augmenting
    % the motif and motifCell
    
    % Adding baseline data on either side, but for analysis we are using
    % only the before BL window for analysis
    % Adding BL times to motifs so it renders in PSTH
    vocalization = [[0 BL_dur]; vocalization]; 
    % Adding BL spikes to motifCell so it renders in PSTH
    motifCell = [cellfun(@(x) x', BL, 'UniformOutput', false) motifCell];
   
    wFs = wavData.fs;
    if isempty(dir(['C:\Users\danpo\Documents\dbh_imgs\*' key '.pdf']))
        
        obj.saveMotifFigs(key, [zeros(BL_buffer * wFs, 1); wavData.wav; zeros(BL_buffer * wFs, 1)],...
            wavData.fs, obj.db(key).amplifier_sampling_rate, motifCell, vocalization, s);
    end
        %}
end

