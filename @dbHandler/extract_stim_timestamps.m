function [timestamps, filecell, adc_sr] = extract_stim_timestamps(obj, S, varargin)

%     Logic:
%     1. Load wav files, text files and audio adc files
%     2. Downsample them to match the recording SR
%     3. Cross correlate signals
%     4. Detect peaks within TTL pulses
%{
%% for re-constructing a text file
for i=1:160
obj.showSpectrogram(norm_adc_audio(over(i):under(i)), S.adc_sr)
input(num2str(i))
cla
end
%}
%{
for i = 1:length(over)
text(over(i), .3, num2str(i))
end
for i = 1:length(under)
text(under(i), -.3, num2str(i))
end
%}

%% 1 set up adc data
disp('adc data')
% adc data
S.board_adc(1,:) = obj.filter_song(S.board_adc(1,:), S.adc_sr);
adc_sr = S.adc_sr;

audioPath = obj.get_audioPath(pwd);

disp('diff')
overunder = [diff(S.board_adc(2,:)) 0]; % add 0 at the end bc this cuts one off

disp('cleaning peaks')
over = find(overunder > 0.95);
under = find(overunder < -0.95);
over = obj.clean_peaks(over);
under = obj.clean_peaks(under);        

% adc audio
disp('adc audio')
adc_audio = S.board_adc(1,:);
norm_adc_audio = adc_audio - mean(adc_audio(1:S.adc_sr));
[b1, a1] = butter(3, 500/S.adc_sr * 2, 'high');
norm_adc_audio =  filter(b1, a1, norm_adc_audio);
[b1, a1] = butter(3, 6e3/S.adc_sr * 2, 'low');
norm_adc_audio =  filter(b1, a1, norm_adc_audio);

% text files for stim id's
if isempty(varargin)
    text_path = dir('*markers.txt');
else
    text_path = dir(varargin{1});
end
fid = fopen(fullfile(text_path(1).folder, text_path(1).name),'r');
filecell = textscan(fid, '%s', 'Delimiter', '\n');
fclose(fid);

if length(over) ~= length(under) || (length(over) ~= length(filecell{1}))
    error('mismatch in number of stim onset and offset times')
end
%% 2 set up audio files
disp('reading audio file')            
wav_files=obj.get_audio(audioPath);

assert(~isempty(wav_files),...
    'You have not added information for this subject to the audioPath property in the dbHandler constructor')
for i=1:length(wav_files)
    cur_wav=wav_files(i).wav;
    cur_fs = wav_files(i).fs;
    cur_wav_resampled = resample(cur_wav, S.adc_sr, cur_fs);
    cur_wav_resampled =  filter(b1, a1, cur_wav_resampled);

    wav_files(i).wav = cur_wav_resampled;
end

%% 3 cross-correlate
%             timestamps = over;
%             hold on;
%             plot(S.board_adc(2,:))
%             plot(S.board_adc(1,:))
%     timestamps = NaN(length(over),1);

for i=1:length(over)
    % extracting the filename pertinent
    curr_wav_path = split(filecell{1}{i},'\');
    curr_wav_name = curr_wav_path{end};
    s = wav_files(strcmp({wav_files.name},curr_wav_name));
%               plot(under(i), 1 ,'r*')
%               plot(under(i)-length(s.data), 1 ,'r*')

    disp(['XCorr Stim: ' int2str(i)])
    % if it's STRF time
    y=[];
    if isempty(s)
        disp(curr_wav_name);
        params = split(curr_wav_name, ' ');
        freq = str2double(params{1});
        amp = str2double(params{2});
        sr = 30e3;
        y = tones(freq, amp, sr);
        % + 1e3 is a buffer because matlab messed up and didn't stop
        % playing before closing the TTL envelope! What the hell!
        [co, lag] = xcorr(norm_adc_audio(over(i):under(i) + 1e3), y); 
    else
        % may need to extend beyond under by two to four seconds
        [co, lag] = xcorr(norm_adc_audio(over(i):under(i)), s.wav);
    end

    [~,I] = max((co));
    lagDiff = lag(I); % is the difference in start of signal between orignial .wav and in TTL envelope
    
    %if ~isempty(s)
    %else
    timestamps(i) = over(i) + (lagDiff); %     under(i) - length(s.wav) / s.fs * adc_sr; % 
    %end
%         Viz
%         figure; hold on; plot(norm_adc_audio(over(i)+lagDiff:under(i)+lagDiff));plot(s.wav/4+.5)
%         figure; hold on; plot(norm_adc_audio(over(i)+lagDiff:under(i)+lagDiff));plot(y/120 +.5)

%         % raw and onoff
    if i == 1
        figure; plot(S.board_adc(2,1:end/3))
        hold on; plot(norm_adc_audio(1:end/3))
    end
    plot(timestamps(i), 0, 'k*')
%                     plot(norm_adc_audio(over(i):under(i)))
%                     plot([zeros(abs(lag(I)),1); y])
%                     plot(lag(I), .2, 'r*')
%                 close all
end
%{
        hold on;
        plot(norm_adc_audio)
        plot(timestamps, ones(length(timestamps),1)*0, 'v');
        hold off;
        close all
%}
end
