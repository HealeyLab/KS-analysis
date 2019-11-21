function [habit_fam, song_fam, pairs, latency_db] = waveform_connector(obj, habit_key, song_key)
%WAVEFORM_CONNECTOR Allows you to take a day of habituation and song and
%connect their respective waveforms
%   Detailed explanation goes here.

%% Get the key families of both inputs
habit_fam = obj.get_key_family(habit_key);
song_fam  = obj.get_key_family(song_key);

%% Generate two figures: one for habit_fam and song_fam
% For both figures, sort the families by channel
    function fam = show_by_channel(obj, fam, above)
        % will sort and plot
        chan_fam = fam';
        for sort_ind = 1:length(fam)
            s = obj.db(fam{sort_ind});
            chan_fam{sort_ind,2} = s.channel;
        end
        [~,I] = sort([chan_fam{:,2}]);
        
        % this value adjusts the axis position for plotting below. 
        % It's either zero (top) or length(fam), which puts it below.
        if above
            row  = 0;
            label = 'habit';
        else
            row = length(fam);
            label = 'song';
        end
        
        % now use these indices to re-arrange the original fam
        for sort_ind = 1:length(fam)
            
            curr_key = chan_fam{I(sort_ind),1};
            chan = chan_fam{I(sort_ind),2};
            
            fam{sort_ind} = curr_key;
            ax = subplot(2,length(fam), sort_ind+row);
            obj.plot_waveform(curr_key, ax);
            % the title is the channel number.
            title(chan);
            
            % in the figure, add data about sort_ind (for later use in this
            % code), and about the tetrode.
            text(ax, 10, 100, ['ind' num2str(sort_ind) ', tet' num2str(floor(chan/4))])
            
            % label the row
            if sort_ind == 1
                ylabel(label);
            end
        end 
    end
figure;
habit_cell = show_by_channel(obj, habit_fam, 0);
song_cell  = show_by_channel(obj, song_fam, 1);

%% Use input to 
done = false;
pairs = {};
while ~done
    user_in = input('top row then bottom row, e.g. "3 3"\n', 's');
    if contains(user_in, 'done')
        done = true;
    else
        user_in = strsplit(user_in, ' ');
        % user_in: first is habit, second is song
        pairs{length(pairs) + 1} = [str2num(user_in{1}) str2num(user_in{2})];
    end
end
%% song data
% DO NOT PUSH THIS TO GITHUB, the original is all I have to fall back on
song_key  = char(song_cell(pairs{1}(1)));
obj.get_song_syllable_activity(song_key, 0); % the one means only that cell, zero means all cells   
input('press enter once youve set the song_syllables','s')
song_sTs = load('C:\Users\danpo\Documents\song.mat'); % song_sTs.song
% Both this function and get_song_syllable_activity use get_key_family.
latency_db = load('C:\Users\danpo\Documents\latency_db.mat');
for i=1:length(pairs)
    %% habit data
    habit_key = char(habit_cell(pairs{i}(1))); % 1 corresponds to habit_cell
    struct_sTs = obj.gen_playback_syllable_PSTHs(habit_key); 
    
    %% find the song_sTs field that best matches song_key
    song_key = char(song_cell(pairs{i}(2))); % 2 corresponds to song
    song_key_rep = strrep(strrep(song_key, ' ', '_'), '&', '_');
    fields = fieldnames(song_sTs.song);
    field = fields{find(strcmp(fieldnames(song_sTs.song), song_key_rep))};
    struct_sTs.song = song_sTs.song.(field);
    
    %% add it to the database
    habit_key = strsplit(habit_key, '_');
    habit_key = strjoin({habit_key{2:end}},'_');
    song_key = strsplit(song_key, '_');
    song_key = strjoin({song_key{2:end}}, '_');
    
    latency_db_key = [habit_key '+' song_key];
    latency_db_key = strrep(latency_db_key, '&', '_'); % ampersands not allowed as field names
    latency_db_key = erase(latency_db_key, 'Kilosort'); % need to shorten it
    latency_db_key = matlab.lang.makeValidName(latency_db_key);
    % initialize as struct so you can call it another field deep later
    latency_db.(latency_db_key) = struct_sTs;
    save('C:\Users\danpo\Documents\latency_db.mat', 'latency_db');    
end
end