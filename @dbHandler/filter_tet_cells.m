function [tet_cells] = filter_tet_cells(obj, tet_cells)
            % removes tet_cells with certain criteria, such as having a low
            % evoked firing rate
            % convert to keys, then to arraylist, then to iterator
            db = obj.db;
%             jList = java.util.ArrayList;jList.add(tet_cells);jList = jList.get(0);it = jList.iterator;
            i = 1;
            while i <= size(tet_cells,1) %for i = 1:length(tet_cells)
                key_i = tet_cells(i,:);
                key = obj.keyhash(key_i{1}, key_i{2}, key_i{3}, key_i{4});
                s = db(key);
                
                if isfield(s, 'stim_timestamps')
                    [on, off, wav_files] = obj.get_stim_on_off(key);

                    si = s.stim_identities{1};
                    usi = unique(si);
                    usi_leg = unique(si);
                    for j = 1:length(usi_leg)
                        edit = usi{j}; edit = strsplit(edit, '\');
                        edit = edit{end}; edit = edit(1:end-4);
                        usi_leg{j} = edit;
                    end

                    fr = 0;
                    % for each class of stim:
                    fr_arr = []; % will overlay response to each of four stims
                    for j = 1:length(usi)
                        stim_inds = find(strcmp(si,usi{j}));
                        % for each stim itself
                        for k = 1:length(stim_inds)
                            sp_ts = s.spike_timestamps;
                            ind = stim_inds(k);
                            evoked = length(intersect(...
                                sp_ts(sp_ts >= on(ind)),...
                                sp_ts(sp_ts < off(ind)))) / ((off(ind)-on(ind)) / s.amplifier_sampling_rate);                      

                            baseline_length_samples = 3 * s.amplifier_sampling_rate;
                            baseline_length_buffer  = 1 * s.amplifier_sampling_rate;

                            baseline = length(intersect(...
                                sp_ts(sp_ts >= (on(ind) - baseline_length_samples)),...
                                sp_ts(sp_ts < on(ind)-baseline_length_buffer))) / (3-1);
                            fr = (evoked - baseline);
                            fr_arr = [fr_arr fr];
                        end
                    end
                    if  mean(fr_arr) < 2 % now in FR!
                        tet_cells(i,:) = [];
                        i = i - 1;
                    end
                end
            i = i + 1;    
            end
        end
        