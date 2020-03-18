function set_song
%%SET_SONG Extracts BOS from an intan recording.
%    This is the function to use to extract BOS at the beginning of an
%    experiment
    
    hs = addcomponents;
    function hs = addcomponents                                             %   [from_left from_bottom width height]
        hs.fig = figure('Visible', 'on', 'Tag', 'fig', 'Units', 'Normalized', 'Position', [0.01 0.06 .95 .8]);
        
        % open files
        hs.open = uicontrol(hs.fig,...
            'Style', 'pushbutton',...
            'Units', 'Normalized',...
            'Position',[0.01, 0.01 0.06 0.04],...
            'String','Open',...
            'UserData', struct('board_adc',[],'Fs',NaN),...
            'Callback', @open,...
            'Tag', 'open');
        
        % same line
        hs.subjectString = uicontrol(hs.fig,...
            'Style', 'pushbutton',...
            'Units', 'Normalized',...
            'Position',[0.05 0.01 0.3 0.04],... % same line
            'Style', 'text',...
            'Units', 'Normalized',...
            'String','Subject name' ); % ("mda", "mdb", ...)'
        
        hs.subjectEdit= uicontrol(hs.fig,...
            'Units', 'Normalized',...
            'Position',[0.25 0.01 0.05  0.04],... % same line
            'Style', 'edit',...
            'String', 'md?',...
            'Tag', 'subject');
        
        % Set the save path
        hs.savePathStr= uicontrol(hs.fig,...
            'Units', 'Normalized',...
            'Style', 'text',...
            'Position',[0.33 0.01 0.09 0.04],... 
            'String','Save path: ');
        hs.savePathEdit = uicontrol(hs.fig,...
            'Units', 'Normalized',...
            'Style', 'edit',...
            'Position',[0.4 0.01 0.3 0.04],...
            'String', 'C:\Users\danpo\Documents\MATLAB\ephysSuite\',...
            'Tag', 'savePath');
        
        % save current selection
        hs.save = uicontrol(hs.fig,...
            'Units', 'Normalized',...
            'Style', 'pushbutton',...
            'Position',[0.9 0.01 0.05  0.04],...
            'String','Save',...
            'Callback', @save,...
            'Tag', 'save');        

        %         [from_left from_bottom width height]
       
        hs.window = axes('Units','Normalized', 'Position', [0.04 0.2 0.95 0.75]);
    end

    function save(~, ~)
        OpenUD = hs.open.UserData;
        Fs = OpenUD.Fs;
        board_adc = OpenUD.board_adc;
        try
            xl = xlim * Fs * 60; % if window needs to convert minutes to seconds
            y = board_adc(end, xl(1):xl(2));

        catch
            xl = xlim * Fs; % if already in seconds,
            y = board_adc(end, xl(1):xl(2));

        end
        y = y-mean(y);
        plot(y)
        audiowrite('BOS.wav', y, Fs);
    end

    function open(hObject, ~)
        [files, origDataPath] = ... % crucial distinction: files vs file (current)
            uigetfile('*.rhd', 'Select an RHD2000 Data File', 'MultiSelect', 'on');

        board_adc = [];
        if ismatrix(files) && ~iscell(files)
            % if is one file
            [~, frequency_parameters, board_adc_data] = read_Intan_RHD2000_file_MML_DJP(fullfile(origDataPath,files),0);
            board_adc = [board_adc board_adc_data];
        else
            % if is two+ files
            filearray = [];
            for i = 1:length(files)
                name = char(files(i));
                path = fullfile(origDataPath, name);
                filearray = [filearray dir(path)];
            end

            [~, idx] = sort({filearray.date});
            files = files(idx);
            
            for i=1:length(files)
                [~, frequency_parameters, board_adc_data] = read_Intan_RHD2000_file_MML_DJP(fullfile(filearray(i).folder,...
                    filearray(i).name),0);
                board_adc = [board_adc board_adc_data];
            end
        end
        
%         S = load(fullfile(origDataPath,'adc_data.mat'));
        hObject.UserData = struct('board_adc', board_adc, 'Fs', frequency_parameters.board_adc_sample_rate); % 
        
        mic = obj.filter_song(board_adc, frequency_parameters.board_adc_sample_rate);
        
        spectrogram(mic, 256, [],[], frequency_parameters.board_adc_sample_rate, 'yaxis')
        colorbar('delete');
    end
end