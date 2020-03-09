function fig_out = gen_song_syllable_PSTHs(obj, key)
% GEN_SONG_SYLLABLE_PSTHS Takes syl_arr as input and outputs psth suubplots
syl_arr = obj.db(key).syl_arr;

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
end