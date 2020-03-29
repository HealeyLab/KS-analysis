function [pbT songT] = extractTables(obj)
%% Initialize
% Pipeline steps: dbh.add, (dbh.waveform_analysis), dbh.waveform_connector
% to connect the habituation and song units. dbh.basicViewer adds the
% timestamps for all the song units for a particular day, and finally
% extractTables puts everything together.

% TODO: subtract 40 ms, simple cellfun I htink
pbT = table('Size', [1 11], 'VariableTypes',...
{'string','double','double',    'double','double','double','double','double','double','double','double'},'VariableNames',...
{'key'   ,'p2p',   'hemisphere','sym',   'Z',     'FR',    'spk',   'lat',   'win',   'BL_FR', 'BL_spk'});
% 1        2        3            4       5        6        7         8        9       10       11

songT = pbT;
%% make tables
% Set intro notes and motifs for pbT, do ths before actually please
% dbh.basicViewer('path\to\md?\BOS.wav')
%%
pairs = readtable('C:\Users\danpo\Documents\MATLAB\DJP_KiloSort\KS-analysis\pairs.xlsx');
for h = 1:height(pairs)
    %% organize data into rectangular motif cell array
    key = pairs.h{h};
    [hs, rasterCell] = obj.gen_psth(key); % rasterCell will be multi line, in seconds
%     axes(hs.sa);
%     drawLines(motifs+2,'y');drawLines(intros+2,'y');
    [introKey, motifKey] = obj.getMotifKey(obj.get_subject_id(key));
    motifs = obj.db(motifKey); intro = obj.db(introKey);
    
    as = obj.get_audio(obj.get_audioPath(key));
    wavData = as(strcmp({as.name}, 'BOS.wav'));
    motifCell = cell(length(rasterCell), size(motifs,1));
    BL = cell(length(rasterCell), 1);BLRange = [0 1.5];
    for r = 1:length(rasterCell)
        % Each row is a different playback. Each column is a motif within
        % that playback
        % [number of distinct motifs, t] = size(motifs)
        for m = 1:size(motifs,1)
            % get timestamps for motifs
            range = motifs(m,:) + 2; % accounting for the 2 seconds added on, I'd like to fix this in the future;
            motifCell{r,m} = obj.sliceTS(rasterCell{r}, range); 
            motifCell{r,m} = motifCell{r,m}'; % needs to be transposed for plotSpikeRaster
        end
        % get BL for pback, will be different for song
        BL{r} = obj.sliceTS(rasterCell{r}, BLRange);
    end
    % finally, call helper
    pbT   = obj.extractTableHelper(pbT, key, motifCell, motifs,...
        wavData.wav, wavData.fs, BL, BLRange);

    %% organize data into nice rectangular array
    key = pairs.s{h};
    [introKey, motifKey] = obj.getMotifKey(key);
    motifs = obj.db(motifKey); intro = obj.db(introKey);
    sTs = obj.db(key).spike_timestamps;
    rasterCell = cell(1, size(motifs,1));
    adc_sr = obj.db(key).adc_sampling_rate;
    for ii = 1:size(motifs,1)
        % only one row, each column is a motif.
        rasterCell{ii} = obj.sliceTS(sTs/adc_sr, motifs(ii,:));
    end
    BL = '????????';
    songT = obj.extractTableHelper(songT, key, rasterCell, motifs,...
        obj.get_microphone(key),adc_sr, BL, BLRange);
end
%% FIX alignments with BOS playbacks if needed

end


