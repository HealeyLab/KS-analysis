function out=genMotifPSTH(obj, sonogram, adc_fs, amp_fs, motifCell, motifs)
%% GENMOTIFPSTH Iteratively runs generatePSTH for each motif        
        for motInd = 1:size(motifs,1)
            mot = motifs(motInd,:) * adc_fs; % converting to samples
            
            % avoiding off-by-one error, converting time zero to index one
            if mot(1) == 0
                mot(1) = 1;
            end
            
            if all(size(motifCell) == [1 1])
                obj.generatePSTH(motifCell, amp_fs,...
                     sonogram(mot(1):mot(2)), adc_fs, size(motifs,1), motInd);
            else
                if mot == [0 0]
                    soundIn = ones(adc_fs,1); % second of silence
                else
                    soundIn = sonogram(mot(1):mot(2));
                end
                
                obj.generatePSTH(motifCell(:,motInd), amp_fs,...
                     soundIn, adc_fs, size(motifs,1), motInd);
            end
        end
        out = gcf;
    end