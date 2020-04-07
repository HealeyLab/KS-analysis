function output = waveform_connector(obj, habit_key, song_key)
%WAVEFORM_CONNECTOR Allows you to take a day of habituation and song and
%connect their respective waveforms
%   Takes the key families of both input arguments and allows you to
%   compare each unit from both recordings.
%%
close all;
%% Get the key families of both inputs
habit_fam = obj.get_key_family(habit_key);
song_fam  = obj.get_key_family(song_key);
polytrode = input('tetrode or stereotrodes?','s');
if contains(polytrode, 'tet')
    polytrode = 4;
elseif contains(polytrode, 'stereo')
    polytrode = 2;
end
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
           
            title(['i:' num2str(sort_ind) '; t:' num2str(floor((chan-1)/polytrode))])
            
            text(ax, 10, 100, [num2str(obj.get_p2p(obj.db(curr_key)))])
            
            % label the row
            if sort_ind == 1
                ylabel(label);
            end
        end 
    end
figure;
habit_cell = show_by_channel(obj, habit_fam, 1);
song_cell  = show_by_channel(obj, song_fam, 0);

%% Use input to connect individual units
done = false;
pairs = {};
while ~done
    user_in = input('Enter the ind of the top row then bottom row, e.g. "3 3"\n', 's');
    if contains(user_in, 'done')
        done = true;
    else
        user_in = strsplit(user_in, ' ');
        % user_in: first is habit, second is song
        pairs{length(pairs) + 1} = [str2num(user_in{1}) str2num(user_in{2})];
    end
end

output = cell(length(pairs),2);
for i=1:length(pairs)
    output{i,1} = habit_cell{pairs{i}(1)};
    output{i,2} = song_cell{pairs{i}(2)};
end

end