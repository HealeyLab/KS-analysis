function hs = generatePSTH(obj, raster_arr,amp_sr, spectData, adc_sr,  varargin)
% range is in samples.
    if size(raster_arr,1) == 1
        ROWS = 2;
    else
        ROWS = 3;
    end
    COLUMNS = 1;
    CURRENT_COLUMN = 1;
    if ~isempty(varargin)
        COLUMNS = varargin{1};
        CURRENT_COLUMN = varargin{2};
    end
        
    start = range(1); stop = range(2);
    hs = struct;
    
    hs.sa = subplot(ROWS, COLUMNS, CURRENT_COLUMN);
    obj.showSpectrogram(spectData, adc_sr);
    
    if ROWS == 3
        subplot(ROWS, COLUMNS, COLUMNS + CURRENT_COLUMN);
    elseif ROWS < 3
        subplot(ROWS, COLUMNS, (ROWS-1) * COLUMNS + CURRENT_COLUMN);
        
    end
    
    [xpoints, ~] = plotSpikeRaster(raster_arr,...
            'PlotTYpe','vertline', 'XLimForCell', [0 length(spectData) / amp_sr]);
        
    hs.ra = subplot(ROWS, COLUMNS, (ROWS-1) * COLUMNS + CURRENT_COLUMN);
    bin = 0.010; % bin size in s
    histogram(hs.ra, xpoints, (0:bin:length(spectData)/amp_sr)); % convert ms to s
    xlim([0, length(spectData)/amp_sr])
end