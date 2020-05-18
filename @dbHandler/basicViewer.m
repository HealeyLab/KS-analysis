function [] = basicViewer(obj, argin, varargin)
%BASICVIEWER Because MATLAB isn't good at OOP, this function extracts
%syllables for both song trials (when you recorded a brain and singing at
%the same time) and for playbacks (for annotating a .wav file to extract
%timestamps
%   To run it for playback, use a .wav absolute path as input and include a
%   second argument with the subject name, to run it for a song trial, just
%   input the entry.
label = 'playback';
if ~isempty(varargin)
    [y, Fs] = audioread(argin);
    sonogram = y;
    adc_sr = Fs;
    argin = varargin{1}; % this will become the key
elseif contains(argin, 'Kilosort')
    % for state of function
    label = 'song';
    %
    entry = obj.db(argin);
    mic = obj.get_microphone(argin); % entry.microphone;
    adc_sr = entry.adc_sampling_rate;
    % Filter
    sonogram = obj.filter_song(mic, adc_sr);
end
%% set up
hs = addcomponents();
axes(hs.sa);

% Midsonogram = medfilt1(sonogram, 5); % Smooth signal
% badIndexes = abs(sonogram-Midsonogram ) > 0.017; % Find bad ones
% sonogram(badIndexes) = Midsonogram(badIndexes); % Repair signal

obj.showSpectrogram(sonogram, adc_sr);
%% declare function variables
intros = [];
motifs = [];


function hs = addcomponents()
    %ADDCOMPONENTS Adds basic components
    %   output is hs, contains handle for parts of gui
    % I have to add an extra layer of helper functions to the
    % callback so I can actually access the hs variable. It's
    % frustrating, but there's nothing I seem to be able to do
    % about it.
    hs=struct;
    hs.fig = figure('Visible', 'on', 'Tag', 'fig', 'Units', 'Normalized', 'Position', [0.01 0.2 .95 .7]);
    hs.sa = axes('Units','Normalized', 'Position', [0.04 0.45 0.95 0.50]);
    hs.ra = axes('Units','Normalized', 'Position', [0.04 0.15 0.95 0.25]);
    hs.btn1 = uicontrol(hs.fig,...
        'Style', 'pushbutton','Units', 'Normalized',...
        'Position',[0.01, 0.02 0.06 0.04],...
        'String','btn1', 'UserData', [],...
        'Callback', @btn1CB, 'Tag', 'btn1');
     hs.btn2 =  uicontrol(hs.fig,...
        'Style', 'pushbutton','Units', 'Normalized',...
        'Position',[0.06, 0.02 0.06 0.04],...
        'String','btn2', 'UserData', [],...
        'Callback', @btn2CB, 'Tag', 'btn2');
    hs.btn3 =  uicontrol(hs.fig,'Style', 'pushbutton','Units', 'Normalized',...
        'Position',[0.17, 0.02 0.06 0.04],...
        'String','sound', 'UserData', [],'Callback', @play, 'Tag', 'open');
    hs.show = uicontrol(hs.fig,...
        'Units', 'Normalized', 'Style', 'pushbutton',...
        'Position',[0.9 0.01 0.06  0.04],...
        'String','show','Callback', @show,...
        'Tag', 'show');  

    function btn1CB(hObject,~)
        intros = [intros; selectWindow('intro')'];
    end

    function btn2CB(hObject,~)
        motifs = [motifs; selectWindow('motif')'];
    end
    
    function play(hObject,~)
        toPlay =sonogram-mean(sonogram);
        axes(hs.sa)
        xl = xlim * adc_sr;
        toPlay = toPlay(xl(1):xl(2));
        sound(toPlay, adc_sr/2);
    end
    
    function show(hObject,~)
        showHelper()
    end
end
   
function showHelper()
    [introKey motifKey]=obj.getMotifKey(argin);
    obj.db(introKey) = intros;
    obj.db([motifKey]) = motifs;
    disp('added')
%       obj.save_db; % do this yourself
end
function x = selectWindow(id)
    [x,~]=ginput(2);
    hs = hs;
    axes(hs.sa);
    line([x(1), x(1)], [0,15e3], 'Color', 'k');
    line([x(2), x(2)], [0,15e3], 'Color', 'k');

    % make label
    axes(hs.sa);
    text(hs.sa, x(1), 12, id, 'FontSize', 18);
end


end


