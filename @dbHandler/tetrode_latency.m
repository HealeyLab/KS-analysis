function latency_map = tetrode_latency(obj)
            % First, separate all keys into 
            keys = obj.db.keys;
            db = obj.db;
            keycell = cell(length(keys), 4 ); % rows, columns
            
            latency_map = containers.Map('KeyType','char','ValueType', 'any');
            
            for i = 1:length(keys)
                % will this present an issue because of unused cell indices
                if isfield(db(keys{i}), 'stim_timestamps')
                    [folder, unit, channel, goodness] = obj.dehash(keys{i});
                    keycell{i,1} = folder; keycell{i,2} = unit;
                    keycell{i,3} = channel; keycell{i,4} = goodness;
                end
            end
            
            keycell = keycell.';
            keycell = reshape(keycell(~cellfun(@isempty,keycell)),4,[])';
            uniquecell = unique(keycell(:,1));
            
            % for each cell, for each recording, for each tetrode...
            for i = 1:length(uniquecell)
                % find cells from same recording
                recording_indices = find(strcmp(uniquecell{i},keycell(:,1)));
                % find cells from same tetrode
                recording_cells = keycell(recording_indices,:); % cells from same recording

                for j = 1:4
                    high = j * 4;
                    low = high - 3;

                    bs = [];
                    ns = [];
                    curr_s = struct('bs', bs, 'ns', ns);


                    % for each tetrode, see if there are shared cells
                    tetrode_cells = recording_cells(...
                        cellfun(@str2num, recording_cells(:,3)) >= low &...
                        cellfun(@str2num, recording_cells(:,3)) <= high,:);

                    for k = 1:size(tetrode_cells,1) % check that this is the right dimension
                        key_i = tetrode_cells(k,:);
                        key_i = obj.keyhash(key_i{1}, key_i{2}, key_i{3}, key_i{4});
                        [on, ~, ~] = obj.get_stim_on_off(key_i);
                        s_i = db(key_i);
                        sp_ts = s_i.spike_timestamps;
                        
                        si = s.stim_identities{1};
                        usi = unique(si);   
                        for l = 1:length(usi)
                            % TODO: add shape for all stimuli
                            all_greater = sp_ts(sp_ts >= on(ind));
                            next_greater = all_greater(1);

                            if s_i.p2p >= 0.43
                                bs = [bs; (next_greater-on(ind))/s.amplifier_sampling_rate];
                            else
                                ns = [ns; (next_greater-on(ind))/s.amplifier_sampling_rate];
                            end

                            curr_s.bs = bs;
                            curr_s.ns = ns;
                            latency_map(key_i) = curr_s;
                        end
                    end
                end
            end
        end
            