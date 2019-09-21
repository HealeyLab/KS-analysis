% same as above function but just for one key
        function fig = cross_correlogram_for_gui(obj,gui_key)
            % First, separate all keys into 
            keys = obj.db.keys;
            keycell = cell(length(keys), 4 ); % rows, columns
            for i = 1:length(keys)
                [folder, unit, channel, goodness] = obj.dehash(keys{i});
                keycell{i,1} = folder;
                keycell{i,2} = unit;
                keycell{i,3} = channel;
                keycell{i,4} = goodness;
            end
            
            [filefolder,~,channel,~] = obj.dehash(gui_key);
            % find cells from same recording
            recording_indices = find(strcmp(filefolder,keycell(:,1)));

            % find cells from same tetrode
            recording_cells = keycell(recording_indices,:); % cells from same recording
            % note: channels can be 1 thru 16, use modulus math to figure
            % out which range a particular one you're dealing with
            
            low = floor(channel / 4) + 1;
            high = low + 3;
            
            tetrode_cells = recording_cells(...
                cellfun(@str2num, recording_cells(:,3)) >= low &...
                cellfun(@str2num, recording_cells(:,3)) <= high,:);

            if size(tetrode_cells,1) > 1
                fig = obj.cross_corr_helper(tetrode_cells);
            end
        end
        