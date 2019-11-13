function fig_out = get_song_activity(obj, key, only)
    % NOTE: to go fast, toggle the fields with the TAB key.
    % add components

    hs = addcomponents;
    show_spectrogram(key);
    syl_arr = [];
    function hs = addcomponents                                             %   [from_left from_bottom width height]        
        
        hs.fig = figure('Visible', 'on', 'Tag', 'fig', 'Units', 'Normalized', 'Position', [0.01 0.2 .95 .7]);
        fig_out = hs.fig;
        hs.A =  uicontrol(hs.fig,...
            'Style', 'pushbutton','Units', 'Normalized',...
            'Position',[0.06, 0.02 0.06 0.04],...
            'String','A', 'UserData', [],...
            'Callback', @make_A, 'Tag', 'open');
        hs.B =  uicontrol(hs.fig,...
            'Style', 'pushbutton','Units', 'Normalized',...
            'Position',[0.11, 0.02 0.06 0.04],...
            'String','B', 'UserData', [],...
            'Callback', @make_B, 'Tag', 'open');
        hs.C =  uicontrol(hs.fig,...
            'Style', 'pushbutton','Units', 'Normalized',...
            'Position',[0.17, 0.02 0.06 0.04],...
            'String','C', 'UserData', [],...
            'Callback', @make_C, 'Tag', 'open');
        hs.D =  uicontrol(hs.fig,...
            'Style', 'pushbutton','Units', 'Normalized',...
            'Position',[0.23, 0.02 0.06 0.04],...
            'String','D', 'UserData', [],...
            'Callback', @make_D, 'Tag', 'open');
        hs.E =  uicontrol(hs.fig,...
            'Style', 'pushbutton','Units', 'Normalized',...
            'Position',[0.29, 0.02 0.06 0.04],...
            'String','E', 'UserData', [],...
            'Callback', @make_E, 'Tag', 'open');
        hs.F =  uicontrol(hs.fig,...
            'Style', 'pushbutton','Units', 'Normalized',...
            'Position',[0.35, 0.02 0.06 0.04],...
            'String','F', 'UserData', [],...
            'Callback', @make_F, 'Tag', 'open');
        hs.G =  uicontrol(hs.fig,...
            'Style', 'pushbutton','Units', 'Normalized',...
            'Position',[0.41, 0.02 0.06 0.04],...
            'String','G', 'UserData', [],...
            'Callback', @make_G, 'Tag', 'open');
        hs.H =  uicontrol(hs.fig,...
            'Style', 'pushbutton','Units', 'Normalized',...
            'Position',[0.47, 0.02 0.06 0.04],...
            'String','H', 'UserData', [],...
            'Callback', @make_H, 'Tag', 'open');
        hs.I =  uicontrol(hs.fig,...
            'Style', 'pushbutton','Units', 'Normalized',...
            'Position',[0.53, 0.02 0.06 0.04],...
            'String','I', 'UserData', [],...
            'Callback', @make_I, 'Tag', 'open');
        hs.J =  uicontrol(hs.fig,...
            'Style', 'pushbutton','Units', 'Normalized',...
            'Position',[0.59, 0.02 0.06 0.04],...
            'String','J', 'UserData', [],...
            'Callback', @make_J, 'Tag', 'open');
        
        %---------------------------------------
        function make_A(hObject,~)
            hObject.UserData = [hObject.UserData make_syllable('A')];
        end
        function make_B(hObject,~)
            hObject.UserData = [hObject.UserData make_syllable('B')];
        end
        function make_C(hObject,~)
            hObject.UserData = [hObject.UserData make_syllable('C')];
        end
        function make_D(hObject,~)
            hObject.UserData = [hObject.UserData make_syllable('D')];
        end
        function make_E(hObject,~)
            hObject.UserData = [hObject.UserData make_syllable('E')];
        end
        function make_F(hObject,~)
            hObject.UserData = [hObject.UserData make_syllable('F')];
        end
        function make_G(hObject,~)
            hObject.UserData = [hObject.UserData make_syllable('G')];
        end
        function make_H(hObject,~)
            hObject.UserData = [hObject.UserData make_syllable('H')];
        end
        function make_I(hObject,~)
            hObject.UserData = [hObject.UserData make_syllable('I')];
        end
        function make_J(hObject,~)
            hObject.UserData = [hObject.UserData make_syllable('J')];
        end
        %---------------------------------
        hs.show = uicontrol(hs.fig,...
            'Units', 'Normalized', 'Style', 'pushbutton',...
            'Position',[0.9 0.01 0.05  0.04],...
            'String','Show','Callback', @show,...
            'Tag', 'show');        

        %         [from_left from_bottom width height]
        hs.align = uicontrol(hs.fig,...
            'Units', 'Normalized', 'Style', 'pushbutton',...
            'Position',[0.8 0.01 0.08  0.04],...
            'String','show syllables','Callback', @align_syllables,...
            'Tag', 'show');        

        hs.sa = axes('Units','Normalized', 'Position', [0.04 0.45 0.95 0.50]);
        hs.ra = axes('Units','Normalized', 'Position', [0.04 0.15 0.95 0.25]);
    end

    function keys = get_keys()
        if only
            keys = key;
        else
            keys = obj.get_key_family(key);
        end
    end

    function align_syllables(~,~)
        % First, organize the syllables according to syllable id
        dict = containers.Map('KeyType','char','ValueType', 'any');
        for syl_ind = 1:length(syl_arr)
            
            curr_syl = syl_arr(syl_ind);
            if ~isKey(dict, curr_syl.id)
                dict(curr_syl.id) = [];                
            end
            dict(curr_syl.id) = [dict(curr_syl.id) curr_syl];
        end
        % Second, if a syllable id has more than one occurrance, take the
        % first occurrance as the basis, and align all the same syllables
        % to it
        % show syllable activity
        % for syllable id in the id,
        dict_keys = dict.keys;
        for i = 1:length(dict_keys)
            sylscells = []; % [rows cells]
            figure('Position', [300, 100, 300, 300]);
            %% for each syllable id:
            curr_syl_id = dict_keys(i); % {'A'}
            curr_syl_id = curr_syl_id{1}; % 'A'
            curr_syls = dict(curr_syl_id); % syllables with the 'A' id
            
            basis_syl = curr_syls(1);
            for j = 1:length(curr_syls)
            %% for each syllable with that syllable id:
            % except for the first, that's the basis
            % currently assuming all syllables are not single...
                lagDiff = 0;
                if j > 1 % if not the basis syllable,
                    curr_syl = curr_syls(j);
                    [co, lag] = xcorr(basis_syl.sonogram, curr_syl.sonogram); 
                    [~,I] = max((co));
                    lagDiff = lag(I); % is the difference in start of signal between orignial .wav and in TTL envelope
                end
                
                subplot(211)
                plot((1:length(curr_syl.sonogram))+lagDiff, curr_syl.sonogram); hold on
                
                % align all syllables to the basis syllable
                syl_cells = curr_syl.cells;
                syl_cells_keys = syl_cells.keys;    
                % redundantly assigned for now
                subplot(212)
                sylscells = [length(curr_syls) length(syl_cells)];
                for k = 1:length(syl_cells_keys)                   
                    %% for each cell of that syllable:
                    % adjust
                    cell_sTs = syl_cells(syl_cells_keys{k}); % get cell sTs
                    cell_sTs = cell_sTs + lagDiff; 
                    
                    % now plot the cell rasters
                    for m = 1:length(cell_sTs)
                        %% for each spike of that cell:
                        % TODO: is this actually plotting different
                        % syllables?
                        % TODO: how are the rasters aligned without
                        % termporal normalization? (ie, subtracting by a
                        % start point
                        y = sylscells(1) * (k-1) + j; % num syllables times cellind norm + cur syl ind
                        rasterRow(cell_sTs(m), y, obj.get_color(obj.db(syl_cells_keys{k}))); 
                    end
                end
            end
            SpectXlim = xlim / curr_syl.adc_sr * curr_syl.amplifier_sr; 
            cla(subplot(211));
            spectrogram(basis_syl.sonogram, 256, [],[], basis_syl.adc_sr, 'yaxis');
            colorbar('delete');
            title(dict_keys{i}); % give it a title
            
            % convert to amplifier sampling rate (for neurons)
            subplot(212)
            xlim(SpectXlim)
            ylim([1 length(curr_syl.cells)*length(curr_syls)+1])
            % shade figure
            hold on;
            num_syls = sylscells(1);
            num_cells = sylscells(2);
            for j = 1:num_cells % total number of regions
                if mod(j, 2) == 1
                    xl = xlim;
                    baseval = (j - 1) * num_syls + 1;
                    h = fill([xl(1) xl(2) xl(2) xl(1)],...
                        [baseval baseval baseval+num_syls baseval+num_syls],...
                        [0.9 0.9 0.9], 'LineStyle', 'None'); 
                    set(h, 'facealpha', 0.5);
                end
            end
        end
    end
    
    function rasterRow(tStamps, i, color)
        line([floor(tStamps),floor(tStamps)], [i, i+1], 'Color', color);
    end
    function show(~, ~)
        cla(hs.ra); % clear the raster
        keys = get_keys(); % get keys
        numUnits = length(keys); % for each cell...
        for i = 1:numUnits
            hold on;
            curr_key = keys{i};
            entry = obj.db(curr_key);

            % adjust xlim
            axes(hs.sa); SpectXlim = xlim * 60 * entry.amplifier_sampling_rate; % convert to amplifier sampling rate (for neurons)
            axes(hs.ra); xlim(SpectXlim)

            % filter timestamps
            sTs = entry.spike_timestamps; % query value
            sTs = sTs(sTs > SpectXlim(1) & sTs < SpectXlim(2)); % only get within this window
            
            color = obj.get_color(entry);
            
            % plot the timestamps
            axes(hs.ra);
            for j = 1:length(sTs)
                rasterRow(sTs(j), i, color);
            end
        end

        % labels
        yticks([1:numUnits]+.5);
        yticklabels([1:numUnits]);
        ylabel('Unit'); xlabel('time (s)');
        axes(hs.ra); prettyAxes(hs.ra, 12);

        axes(hs.sa); xlabel('')
        prettyAxes(hs.sa, 12);
    end
    function prettyAxes(axis, fontSize)
        axis.YAxis.FontSize = fontSize; axis.XAxis.FontSize = fontSize;
    end
    function show_spectrogram(key)
        % remove all keys not in key family
        entry = obj.db(key);
        mic = obj.get_microphone(key); % entry.microphone;
        adc_sr = entry.adc_sampling_rate;

        % Filter
        fcutlow = 300;
        fcuthigh = 12e3;
        order = 3;
        [b,a]=butter(order,[fcutlow,fcuthigh]/(adc_sr/2),'bandpass');
        mic = filter(b, a, mic);
        mic = fliplr(mic);
        mic = filter(b, a, mic);
        mic = fliplr(mic);

        axes(hs.sa);
        spectrogram( mic(end,:), 256, [],[], adc_sr, 'yaxis')
        colorbar('delete');
    end
%  --------------------------------------
    function curr_syl = make_syllable(id)
        [x,~]=ginput(2);
        axes(hs.sa);
        line([x(1), x(1)], [0,15], 'Color', 'k');
        line([x(2), x(2)], [0,15], 'Color', 'k');
        % make label
        xmax = xlim; xmax = max(xmax);
        ymax = ylim; ymax = max(ymax);
        axes(hs.sa);
        text(hs.sa, x(1), 12, id, 'FontSize', 18);
        
        % get spiketrains
        i=0;
        keys = get_keys();
        curr_syl = syllable(obj, x, keys, id);
        syl_arr = [syl_arr curr_syl];
    end
end