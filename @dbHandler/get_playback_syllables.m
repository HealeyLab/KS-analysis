function fig_out = get_playback_syllables(obj)
%GET_SUBJECT_SYLLABLES Takes syllables from the song recordings from
%playback experiments
%   For playback experiments, we need to be able to do syllable-by-syllable
%   analyses. To do that, I'm going to save the syllables for each
%   subject's calls in the database under a custom key (just subject name)
%   

% db key-value pair, where the valye is a struct with fields CON, BOS,...,
% with values for those fields of syllable objects.

subject = 'mda';
path = 'C:\Users\danpo\Documents\MATLAB\ephysSuite\zf son mda';
hs = addcomponents;
syl_arr = [];

syl_struct = struct;
audio_s = obj.get_audio(path);
for i = 1:length(audio_s)
    curr_s = audio_s(i);
    % add to struct
    fieldname = split(curr_s.name, '.');
    fieldname = fieldname{1};
    
    if ~contains(curr_s.name, 'whitenoise')
        curr_s = audio_s(i);
        show_spectrogram(curr_s.wav, curr_s.fs, hs.sa);
        % go to next
        inp = '';
        while ~strcmp(inp, 'continue')
            inp = input('enter continue or play\n', 's');
            if strcmp(inp, 'play')
                sound(curr_s.wav, curr_s.fs/2); % half Shpeed!
            elseif strcmp(inp, 'reset')
                cla(hs.sa)
                syl_arr = [];
                show_spectrogram(curr_s.wav, curr_s.fs, hs.sa);
            end
        end
        syl_struct.(fieldname) = syl_arr;
    end

    % reset syl_arr to empty
    syl_arr = [];
end

subject_key = [subject '_playback_syllables'];
obj.db(subject_key) = syl_struct;

disp('all done')

function show_spectrogram(son, fs, ax)
    % Filter
    fcutlow = 300; fcuthigh = 12e3;
    order = 3;
    [b,a]=butter(order,[fcutlow,fcuthigh]/(fs/2),'bandpass');
    son = filter(b, a, son);son = fliplr(son);
    son = filter(b, a, son);son = fliplr(son);

    axes(ax);
    spectrogram(son(:,end), 512, [],[], fs, 'yaxis')
    colorbar('delete');
end
% assemble the gui
function hs = addcomponents
    hs.fig = figure('Visible', 'on', 'Tag', 'fig', 'Units', 'Normalized', 'Position', [0.01 0.2 .95 .7]);
    
    hs.Intro = uicontrol(hs.fig,...
        'Style', 'pushbutton','Units', 'Normalized','Position',[0.01, 0.02 0.06 0.04],...
        'String','Intro', 'UserData', [], 'Callback', @make_Intro, 'Tag', 'open');
    hs.A =  uicontrol(hs.fig,...
        'Style', 'pushbutton','Units', 'Normalized','Position',[0.06, 0.02 0.06 0.04],...
        'String','A', 'UserData', [],'Callback', @make_A, 'Tag', 'open');
    hs.B =  uicontrol(hs.fig, 'Style', 'pushbutton','Units', 'Normalized',...
        'Position',[0.11, 0.02 0.06 0.04],'String','B', 'UserData', [],...
        'Callback', @make_B, 'Tag', 'open');
    hs.C =  uicontrol(hs.fig,'Style', 'pushbutton','Units', 'Normalized',...
        'Position',[0.17, 0.02 0.06 0.04],...
        'String','C', 'UserData', [],'Callback', @make_C, 'Tag', 'open');
    hs.D =  uicontrol(hs.fig,'Style', 'pushbutton','Units', 'Normalized',...
        'Position',[0.23, 0.02 0.06 0.04], 'String','D', 'UserData', [],...
        'Callback', @make_D, 'Tag', 'open');
    hs.E =  uicontrol(hs.fig,'Style', 'pushbutton','Units', 'Normalized',...
        'Position',[0.29, 0.02 0.06 0.04],'String','E', 'UserData', [],...
        'Callback', @make_E, 'Tag', 'open');
    hs.F =  uicontrol(hs.fig,...
        'Style', 'pushbutton','Units', 'Normalized','Position',[0.35, 0.02 0.06 0.04],...
        'String','F', 'UserData', [],'Callback', @make_F, 'Tag', 'open');
    hs.G =  uicontrol(hs.fig, 'Style', 'pushbutton','Units', 'Normalized',...
        'Position',[0.41, 0.02 0.06 0.04],'String','G', 'UserData', [],...
        'Callback', @make_G, 'Tag', 'open');
    hs.H =  uicontrol(hs.fig,...
        'Style', 'pushbutton','Units', 'Normalized','Position',[0.47, 0.02 0.06 0.04],...
        'String','H', 'UserData', [],'Callback', @make_H, 'Tag', 'open');
    hs.I =  uicontrol(hs.fig,'Style', 'pushbutton','Units', 'Normalized',...
        'Position',[0.53, 0.02 0.06 0.04],'String','I', 'UserData', [],...
        'Callback', @make_I, 'Tag', 'open');
    hs.J =  uicontrol(hs.fig,'Style', 'pushbutton','Units', 'Normalized',...
        'Position',[0.59, 0.02 0.06 0.04],'String','J', 'UserData', [],...
        'Callback', @make_J, 'Tag', 'open');
    hs.backspace = uicontrol(hs.fig,'Units', 'Normalized', 'Style', 'pushbutton',...
        'Position',[0.8 0.06 0.08  0.04],'String','backspace','Callback', @backspace,...
        'Tag', 'show');
    %--------------------------------------------------------
    hs.sa = axes('Units','Normalized', 'Position', [0.04 0.15 0.95 0.75]);

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
    %-----------------------------------------------
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
        curr_syl = syllable(obj, x, {}, id);
        syl_arr = [syl_arr curr_syl];
    end

    function backspace(~,~)
        syl_arr(length(syl_arr)) = [];
    end
end
end