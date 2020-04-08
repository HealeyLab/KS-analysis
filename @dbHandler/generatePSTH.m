function hs = generatePSTH(obj, raster_arr,amp_sr, spectData, adc_sr,  varargin)
%% GENERATEPSTH Does exactly what the function says
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
        
    %% Generate spectrogram
    start = range(1); stop = range(2);
    hs = struct;
    
    hs.sa = subplot(ROWS, COLUMNS, CURRENT_COLUMN);
    obj.showSpectrogram(spectData, adc_sr);
    set(gca,'xticklabel',[])
    set(gca,'yticklabel',[])
    set(gca,'xlabel',[])
    set(gca,'ylabel',[])
    
    %% generate raster
    if ROWS == 3
        subplot(ROWS, COLUMNS, COLUMNS + CURRENT_COLUMN);
    elseif ROWS < 3
        subplot(ROWS, COLUMNS, (ROWS-1) * COLUMNS + CURRENT_COLUMN);
    end
    % When it gives you a transpose error, it's because one of the inner
    % arrays are not tansposed right. It's not the specific cell array,
    % it's one (or all) of the comprising arrays.
    [xpoints, ~] = plotSpikeRaster(raster_arr,...
            'PlotTYpe','vertline', 'XLimForCell', [0 length(spectData) / adc_sr]);
    xpoints = xpoints(1:3:end); % every third value, for some reason it's redundant.
    set(gca,'xticklabel',[])
    set(gca,'yticklabel',[])
    
    %% generate histogram    
    hs.ra = subplot(ROWS, COLUMNS, (ROWS-1) * COLUMNS + CURRENT_COLUMN);
    bin = 0.010; % bin size in s
    duration = length(spectData)/adc_sr;
    count = histcounts(xpoints,floor(duration/bin)); % duration/bin is the number of bins
    bar(count,'k')
    set(gca,'xticklabel',[])
    
    xlabel([num2str(round(duration,3)) 's'])
    
    if CURRENT_COLUMN == 1
        ylabel('FR (Hz)')
    end
    %%
    function drawLines(arr, c)
            hold on;

        for stimInd = 1:size(arr,1)
            line([arr(stimInd,1), arr(stimInd,1)], [0,15e3], 'Color', c);
            line([arr(stimInd,2), arr(stimInd,2)], [0,15e3], 'Color', c);
        end
    end
end