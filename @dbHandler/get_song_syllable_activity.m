function get_song_syllable_activity(obj, key, only)
%%GET_SONG_SYLLABLE_ACTIVITY For directed song recordings, this script is
%%for taking the timestamps of the syllables in the recorded song.
% Saves data in a separate database located in Documents/song.mat
% To go fast, toggle the fields with the TAB key.
%% add components

    hs = addcomponents;
    show_spectrogram(key);
    syl_arr = [];
    % FOR OUTPUT
    song = struct;
    % FOR OUTPUT
    function hs = addcomponents    %   [from_left from_bottom width height]        
        
        hs.fig = figure('Visible', 'on', 'Tag', 'fig', 'Units', 'Normalized', 'Position', [0.01 0.2 .95 .7]);
        fig_out = hs.fig;
        
        hs.Intro = uicontrol(hs.fig,...
            'Style', 'pushbutton','Units', 'Normalized',...
            'Position',[0.01, 0.02 0.06 0.04],...
            'String','Intro', 'UserData', [],...
            'Callback', @make_Intro, 'Tag', 'open');
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
        
        function make_Intro(hObject,~)
            hObject.UserData = [hObject.UserData make_syllable('Intro')];
        end
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
            'Position',[0.9 0.01 0.06  0.04],...
            'String','show whole','Callback', @show,...
            'Tag', 'show');        

        %         [from_left from_bottom width height]
        hs.align = uicontrol(hs.fig,...
            'Units', 'Normalized', 'Style', 'pushbutton',...
            'Position',[0.8 0.01 0.08  0.04],...
            'String','show syllables','Callback', @align_syllables,...
            'Tag', 'show');
        
        hs.backspace = uicontrol(hs.fig,...
            'Units', 'Normalized', 'Style', 'pushbutton',...
            'Position',[0.8 0.06 0.08  0.04],...
            'String','backspace','Callback', @backspace,...
            'Tag', 'show');
        
        hs.sa = axes('Units','Normalized', 'Position', [0.04 0.45 0.95 0.50]);
        hs.ra = axes('Units','Normalized', 'Position', [0.04 0.15 0.95 0.25]);
    end
    function keys = get_keys()
        if only
            keys = {key};
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
        % for syllable id in the id,
        dict_keys = dict.keys;
        figure;
        for i = 1:length(dict_keys)            
            sylscells = []; % [rows cells]
%             figure('Position', [300, 100, 300, 300]);
            %% for each syllable id:
            curr_syl_id = dict_keys(i); % {'A'}
            curr_syl_id = curr_syl_id{1}; % 'A'
            curr_syls = dict(curr_syl_id); % syllables with the 'A' id
            
            basis_syl = curr_syls(1);
            for j = 1:length(curr_syls)
            %% for each syllable with that syllable id:
            % except for the first, that's the basis
                lagDiff = 0;
                curr_syl = curr_syls(j);

                if j > 1 % if not the basis syllable,
                    [co, lag] = xcorr(basis_syl.sonogram, curr_syl.sonogram); 
                    [~,I] = max((co));
                    lagDiff = lag(I); % is the difference in start of signal between orignial .wav and in TTL envelope
                end
                
%                 subplot(...)
%                 plot((1:length(curr_syl.sonogram))+lagDiff, curr_syl.sonogram); hold on
                
                % align all syllables to the basis syllable
                all_cells = curr_syl.cells;
                subplot(2, length(dict), length(dict)+i);
                sylscells = [length(curr_syls) length(all_cells)];
                
                % FOR OUTPUT
                if j == 1 % if this is the basis syllable
                    % because if it's the first iteration, then hasn't been
                    % initialized. Needs to be initialized to something to
                    % assign it a value. Weird struct rule.
                    
                    for k = 1:length(all_cells)
                        if i == 1
                            song.(strrep(strrep(curr_syl.cell_keys{k}, ' ', '_'), '&', '_')) = struct;
                            song.(strrep(strrep(curr_syl.cell_keys{k}, ' ', '_'), '&', '_'))...
                                .(curr_syl_id) = struct;
                        end
                       
                        % this is ugly, I'm replacing the ampersands and spaces
                        % with underscores to initialize all fields as structs.
                        % Love is war.
                        song.(strrep(strrep(curr_syl.cell_keys{k}, ' ', '_'), '&', '_'))...
                            .(curr_syl_id).sTs = {};
                        song.(strrep(strrep(curr_syl.cell_keys{k}, ' ', '_'), '&', '_'))...
                            .(curr_syl_id).window_s = {};
                    end
                end
                % FOR OUTPUT
                
                for k = 1:length(all_cells)                   
                    %% for each cell of that syllable:
                    % adjust
                    cell_sTs = all_cells{k}; % get cell sTs
                    cell_sTs = cell_sTs + lagDiff; 
                    cell_sTs = cell_sTs;
                    
                    % FOR OUTPUT
                    song.(strrep(strrep(curr_syl.cell_keys{k}, ' ', '_'), '&', '_'))...
                        .(curr_syl_id).sTs...
                        = [song.(strrep(strrep(curr_syl.cell_keys{k}, ' ', '_'), '&', '_')).(curr_syl_id).sTs; cell_sTs];
                    song.(strrep(strrep(curr_syl.cell_keys{k}, ' ', '_'), '&', '_'))...
                        .(curr_syl_id).window_s ...
                        = [song.(strrep(strrep(curr_syl.cell_keys{k}, ' ', '_'), '&', '_')).(curr_syl_id).window_s; curr_syl.window_s];
                    
                    % FOR OUTPUT
                    
                    % now plot the cell rasters
                    for m = 1:length(cell_sTs)
                        %% for each spike of that cell:
                        y = length(curr_syls) * (k-1) + j; % num syllables times cellind norm + cur syl ind
                        rasterRow(cell_sTs(m), y, obj.get_color(obj.db(curr_syl.cell_keys{k}))); 
                    end
                end
            end
            SpectXlim = xlim / curr_syl.adc_sr * curr_syl.amplifier_sr; 
%             cla(subplot(2,length(dict),i));
            subplot(2,length(dict),i)
            spectrogram(basis_syl.sonogram, 256, [],[], basis_syl.adc_sr, 'yaxis');
            colorbar('delete');
            title(dict_keys{i}); % give it a title
            
            % convert to amplifier sampling rate (for neurons)
            subplot(2,length(dict),length(dict) + i)
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
        %% save syl_arr to temp file
        save('C:\Users\danpo\Documents\song.mat', 'song')
        % add syl_arr to obj so you can gen psths later!
        s = obj.db(key);
        s.syl_arr = syl_arr;
        obj.db(key) = s;
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
            % 
            sTs = entry.spike_timestamps; % query value
            % bring back in time, align with the audio by subtracting 40 ms
            sTs = sTs - 0.040 * entry.adc_sampling_rate; 
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
        mic = obj.filter_song(mic, adc_sr);

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
        keys = get_keys();
        curr_syl = syllable(obj, x, keys, id);
        syl_arr = [syl_arr curr_syl];
    end
    function backspace(~,~)
        syl_arr(length(syl_arr)) = [];
    end
end