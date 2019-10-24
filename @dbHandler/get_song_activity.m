function fig_out = get_song_activity(obj, key, only)
    % NOTE: to go fast, toggle the fields with the TAB key.
    % add components

    hs = addcomponents;
    function hs = addcomponents                                             %   [from_left from_bottom width height]
        hs.fig = figure('Visible', 'on', 'Tag', 'fig', 'Units', 'Normalized', 'Position', [0.01 0.2 .95 .7]);
        fig_out = hs.fig;
        hs.open = uicontrol(hs.fig,...
            'Style', 'pushbutton','Units', 'Normalized',...
            'Position',[0.01, 0.01 0.06 0.04],...
            'String','Open',...
            'UserData', struct('board_adc',[],'Fs',NaN, 'key',''),...
            'Callback', @open, 'Tag', 'open');

        % same line
        hs.subjectString = uicontrol(hs.fig,...
            'Style', 'pushbutton', 'Units', 'Normalized',...
            'Position',[0.05 0.01 0.3 0.04],... % same line
            'Style', 'text', 'Units', 'Normalized',...
            'String','Subject name' ); % ("mda", "mdb", ...)'

        hs.subjectEdit= uicontrol(hs.fig,...
            'Units', 'Normalized',...
            'Position',[0.25 0.01 0.05  0.04],... % same line
            'Style', 'edit',...
            'String', 'md?', 'Tag', 'subject');

        % Set the save path
        hs.savePathStr= uicontrol(hs.fig,...
            'Units', 'Normalized', 'Style', 'text',...
            'Position',[0.33 0.01 0.09 0.04],... 
            'String','Save path: ');
        hs.savePathEdit = uicontrol(hs.fig,...
            'Units', 'Normalized', 'Style', 'edit',...
            'Position',[0.4 0.01 0.3 0.04],...
            'String', 'C:\Users\danpo\Documents\MATLAB\ephysSuite\',...
            'Tag', 'savePath');

        % save current selection
        hs.save = uicontrol(hs.fig,...
            'Units', 'Normalized', 'Style', 'pushbutton',...
            'Position',[0.9 0.01 0.05  0.04],...
            'String','Show','Callback', @show,...
            'Tag', 'show');        

        %         [from_left from_bottom width height]

        hs.sa = axes('Units','Normalized', 'Position', [0.04 0.45 0.95 0.50]);
        hs.ra = axes('Units','Normalized', 'Position', [0.04 0.15 0.95 0.25]);
    end

    function show(~, ~)
        % shows on multiple lines
        function rasterRow(tStamps, i, color)
            axes(hs.ra);
            line([floor(tStamps),floor(tStamps)], [i, i+1], 'Color', color);
        end

        % clear the raster
        cla(hs.ra);

        OpenUD = hs.open.UserData;
        keys = OpenUD.keys;

        % for each cell
        numUnits = length(OpenUD);
        for i = 1:numUnits
            hold on;
            curr_key = OpenUD(i).keys;
            entry = obj.db(curr_key);

            % get the time limits in samples
            aFs = entry.amplifier_sampling_rate;

            % adjust xlim
            axes(hs.sa); SpectXlim = xlim * 60 * aFs; % right into amplifier sampling rate
            axes(hs.ra); xlim(SpectXlim)

            % filter timestamps
            sTs = entry.spike_timestamps; % query value
            sTs = sTs(sTs > SpectXlim(1) & sTs < SpectXlim(2)); % only get within this window
            
            p2p = obj.get_p2p(entry);
            if p2p >= .43
                color = obj.BB;
            else
                color = obj.NN;
            end
            
            % plot the timestamps
            for j = 1:length(sTs)
                rasterRow(sTs(j), i, color);
            end
        end

        % labels
        yticks([1:numUnits]+.5)
        yticklabels([1:numUnits])
        ylabel('Unit')
        xlabel('time (s)');
        axes(hs.ra)
        prettyAxes(hs.ra, 12);

        axes(hs.sa);
        xlabel('')
        prettyAxes(hs.sa, 12);

    end
    function prettyAxes(axis, fontSize)
        axis.YAxis.FontSize = fontSize;
        axis.XAxis.FontSize = fontSize;
    end

    
    function open(hObject, ~)
        % remove all keys not in key family
        keyOrig = key;
        if only
            keys = key;
        else
            keys = obj.get_key_family(keyOrig);
        end
        entry = obj.db(keyOrig);
        mic = entry.microphone;
        adc_sr = entry.adc_sampling_rate;
        hObject.UserData = struct(...
            'board_adc', mic,...
            'Fs', adc_sr,...
            'keys', keys); % frequency_parameters.board_adc_sample_rate

        axes(hs.sa);
        spectrogram( mic(end,:), 256, [],[], adc_sr, 'yaxis')
        colorbar('delete');
    end
end