function all_sTs = gen_playback_syllable_PSTHs(obj, key)
%% GEN_PLAYBACK_SYLLABLE_PSYTHs
% The output of this function is a heirarchichal database struct whos fields are named after each
% song type. The value for each field is the cell_sTs for 

% key = 'mde 10 13 19 habituation_191013_173207_Kilosort&152&13&good'; % obj.get_keys(subject, '_syllables'); % using & as a pattern excludes the auxiliary entries
obj_keys = {key}; % obj.get_key_family(seed_key);
num_cells = length(obj_keys); % always 1
[subject,~] = obj.get_subject_id(key);

%%
syllable_key = [subject '_playback_syllables'];
subject_songs = obj.db(syllable_key);
fn = fieldnames(subject_songs);

[on, off, ~] = obj.get_stim_on_off(key);
stim_identities = obj.db(key).stim_identities{1};

% FOR OUTPUT
all_sTs = struct;
% FOR OUTPUT
%% for each type of song (BOS, CON, BOS-REV)...
for song_index = 1:length(fn)
    % FOR OUTPUT
    all_sTs.(fn{song_index}) = struct;
    % FOR OUTPUT
    
     [cur_wav, fs] = audioread(fullfile(...
         ['C:\Users\danpo\Documents\MATLAB\ephysSuite\zf son ' subject],...
         [fn{song_index} '.wav']));
    song = subject_songs.(fn{song_index});
    num_syllables = length(song);
    syllable_dict = containers.Map('KeyType','char','ValueType', 'any');
    % build dict
    for i = 1:num_syllables

        curr_syl = song(i);
        if ~isKey(syllable_dict, curr_syl.id)
            syllable_dict(curr_syl.id) = [];
        end
        syllable_dict(curr_syl.id) = [syllable_dict(curr_syl.id) curr_syl];
    end
    % Second, if a syllable id has more than one occurrance, take the
    % first occurrance as the basis, and align all the same syllables
    % for syllable id in the id,
    dict_keys = syllable_dict.keys;
    figure('Name',fn{song_index});
    title(strrep(fn{song_index},'_', ' '));
    for i = 1:length(dict_keys)
        
        %% for each syllable id:
        curr_syl_id = dict_keys(i); % {'A'}
        curr_syl_id = curr_syl_id{1}; % 'A'
        
        syls_in_curr_type_of_song = syllable_dict(curr_syl_id); % syllables with the 'A' id

        basis_syl = syls_in_curr_type_of_song(1);
        basis_syl_sonogram = get_sonogram(cur_wav, fs, basis_syl);
        
        % FOR OUTPUT
        all_sTs.(fn{song_index}).(curr_syl_id) = {};
        % FOR OUTPUT
        for j = 1:length(syls_in_curr_type_of_song)
        %% for each syllable with that syllable id:
        % except for the first, that's the basis
            lagDiff = 0;
            curr_syl = syls_in_curr_type_of_song(j);
            curr_syl_sonogram = get_sonogram(cur_wav, fs, curr_syl);
            if j > 1 % if not the basis syllable,
                [co, lag] = xcorr(basis_syl_sonogram, curr_syl_sonogram); 
                [~,I] = max((co));
                lagDiff = lag(I); % is the difference in start of signal between orignial .wav and in TTL envelope
            end
            
            subplot(2, length(syllable_dict), i)
            plot((1:length(curr_syl_sonogram))+lagDiff, curr_syl_sonogram); hold on

            % align all syllable-evoked activity to the basis syllable
            subplot(2, length(syllable_dict), length(syllable_dict)+i);
            for k = 1:num_cells % aka length(obj_keys)           
                %% for each cell of that syllable:
                % % to test that it works:
                % plot(board_adc(1,:))
                % hold on
                % for ind = 1:length(song)
                %     plot(song(ind).window_s(1)*30000+on_(1), .1, 'r*')
                % end
                total_sTs = obj.db(obj_keys{k}).spike_timestamps;
                % bring back in time to align with sound
                total_sTs = total_sTs - 0.040 * obj.db(obj_keys{k}).amplifier_sampling_rate;
                % timestamp of each stimulus of the current song type
                on_   = on(contains(stim_identities,[fn{song_index} '.wav']));
                off_ = off(contains(stim_identities,[fn{song_index} '.wav']));
                % make a cell each stimulus presentation and
                % subtract the onset of each stimulus from the
                % corresponding cell array entry then, take spikes from
                % that syllable only
                curr_sTs = cell(length(on_),1);
                for m = 1:length(on_)
                    % Each stimulus presentation
                    curr_sTs{m} = total_sTs(...
                          total_sTs >  on_(m)...
                        & total_sTs < off_(m));
                    curr_sTs{m} = curr_sTs{m} - on_(m);
                    % NOTE: BECAUSE SONGS ARE LESS THAN A MINUTE, THE XLIM FOR THE SPECTROGRAM
                    % IS IN SECONDS AND NOT MINUTES, SO ADJUST ACCORDINGLY
                    % for each song rendition,
                    
                    % Each syllable presented
                    amp = sort(curr_syl.window_s * obj.db(obj_keys{k}).amplifier_sampling_rate); % convert to samples. not multiplying by 60 as I have elsewhere
                    cs = curr_sTs{m};
                    curr_sTs{m} = cs(...
                          curr_sTs{m} > amp(1)...
                        & curr_sTs{m} < amp(2));
                    curr_sTs{m} = curr_sTs{m} - amp(1);
                    curr_sTs{m} = curr_sTs{m} + lagDiff;
                end
                
                % FOR OUTPUT
                all_sTs.(fn{song_index}).(curr_syl_id) = [all_sTs.(fn{song_index}).(curr_syl_id); curr_sTs];
                % FOR OUTPUT
                
                % now plot the cell rasters
                for row_ind = 1:length(curr_sTs)
                    row = curr_sTs{row_ind};
                    for sTs_ind = 1:length(row)
                        %% for each spike of that cell:
                        % j is the nth syllable presentation
                        % length(On_) is the total number of syllable
                        % presentations
                        y = (j-1) * length(on_) + (row_ind);
                        rasterRow(row(sTs_ind), y, obj.get_color(obj.db(obj_keys{k}))); 
                    end
                end
            end
        end
        SpectXlim = xlim / fs * obj.db(obj_keys{k}).amplifier_sampling_rate; 
        cla(subplot(2,length(syllable_dict),i));
        spectrogram(basis_syl_sonogram, 256, [],[], fs, 'yaxis');
        colorbar('delete');
        title(dict_keys{i}); % give it a title

        % convert to amplifier sampling rate (for neurons)
        subplot(2,length(syllable_dict),length(syllable_dict) + i)
        xlim(SpectXlim)
        ylim([1 length(on_)*length(syls_in_curr_type_of_song)+1])
        % shade figure
        hold on;
        for j = 1:length(syls_in_curr_type_of_song) % total number of regions
            if mod(j, 2) == 1
                xl = xlim;
                baseval = (j - 1) * length(on_) + 1;
                h = fill([xl(1) xl(2) xl(2) xl(1)],...
                    [baseval baseval baseval+length(on_) baseval+length(on_)],...
                    [0.9 0.9 0.9], 'LineStyle', 'None'); 
                set(h, 'facealpha', 0.5);
            end
        end
    end
end
    function sonogram = get_sonogram(wav, fs, syl)
        window_s = floor(sort(syl.window_s) * fs);
        sonogram = wav(window_s(1):window_s(2));
    end
    function rasterRow(tStamps, i, color)
        line([floor(tStamps),floor(tStamps)], [i, i+1], 'Color', color,'LineWidth',2);
    end

%%
% song presentation
%              1             |              2             |   |  50
% syllables                  |                            |   |
% A1 B1 C1 A2 B2 C2 A3 B3 C3 | A1 B1 C1 A2 B2 C2 A3 B3 C3 |...| A1 B1 C1 A2 B2 C2 A3 B3 C3
%____________________________|____________________________|___|_
%
% raster for A:
%|------------
%|   .     . 
%|  . . . ... .    A3
%|  . .. .  .
%|......
%|------------
%|  .     . 
%|  . . . ... .    A2
%|  . .. .  .
%|......
%|------------
%|. . ... . 
%|.    . .. .      A1 
%|_____________
end