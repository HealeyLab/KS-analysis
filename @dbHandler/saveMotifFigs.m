function saveMotifFigs(obj,key, sonogram, adc_sr, amplifier_sr, motifCell, motifs, varargin)
%                          1    2         3       4             5          6       7 
%ADDMOTIFFIGS Generates and saves raw data figures, then closes them
%{   
    Inputs
        -key: dbHandler key
        - sonogram: time series sound signal that you want used in PSTH
        - adc_sr: sampling rate for sonogram signal
        - amplfier_sr: sampling rate of spikes (ie, conversion factor for
        spikes from seconds to samples)
        - motifCell: M x 1 cell of 1 x N arrays of spike timestamps
        - motifs: M x 2 array of onsets and offsets of motifs
        - varargin: may contain spike timestamps for plotting the 
%} 

key_text = replace(key, '_', ' ');
figure;
obj.genMotifPSTH(sonogram, adc_sr, amplifier_sr, motifCell,motifs); 
% motif PSTH
% addToPDF(gcf, [getSecond key]);

% Playback
if isstruct(varargin{1}) 
    s = varargin{1};
    obj.gen_psth(s, 1);
% Song trial
elseif isnumeric(varargin{1})
    sTs = varargin{1};
    range = [min(min(motifs)) max(max(motifs))]; % in seconds

    obj.generatePSTH({obj.sliceTS(sTs/amplifier_sr, range)},...
        amplifier_sr, sonogram(range(1)*adc_sr:range(2)*adc_sr), adc_sr);
else
    error('difficulty handling varargin')
end    

% whole PSTH
addToPDF(gcf, [getSecond key_text]); 

%% helpers
    function addToPDF(cf, text)
        output = ['C:\Users\danpo\Documents\dbh_imgs\' text '.png'];
        mtit(strrep(text, '_', ' '), 'fontsize', 20)
%         text = strsplit(text,'Kilo');
        set(cf, 'Position', get(0, 'Screensize'));
        export_fig (output)
%         append_pdfs('C:\Users\danpo\Documents\MATLAB\test.pdf', 'C:\Users\danpo\Documents\MATLAB\temp.pdf')
        close(cf); % making room for the next thing
    end

    function secs = getSecond()
        secs=datetime('now');

        secs = datestr(secs);
        secs = strsplit(secs, ':');
        secs  = secs{end};
    end
end

