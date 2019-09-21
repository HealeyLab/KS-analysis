
        function [BBavg, NNavg, BNavg, NBavg] = cross_correlograms(obj)
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
            BBavg = []; NNavg = []; BNavg = []; NBavg = [];
            uniquecell = unique(keycell(:,1));
            
            % find cells from same recording
            for i = 1:length(uniquecell)
                recording_indices = find(strcmp(uniquecell{i},keycell(:,1)));
                
                % find cells from same tetrode
                recording_cells = keycell(recording_indices,:); % cells from same recording
                % note: channels can be 1 thru 16
                % for each tetrode, see if there are shared cells
                for j = 1:4
                    high = j * 4;
                    low = high - 3;

                    tetrode_cells = recording_cells(...
                        cellfun(@str2num, recording_cells(:,3)) >= low &...
                        cellfun(@str2num, recording_cells(:,3)) <= high,:);
                    
                    if size(tetrode_cells,1) > 1
                        xcorr_filename = ['C:\Users\danpo\Documents\dbh_imgs\'...
                            uniquecell{i} '_tetrode' num2str(floor(low/4) + 1) '_xcorr_filtered.png']; % previously saved as jpg, but now I need high quality stuff 
                        
                        [fig, tosave,...
                        BBavg,NNavg,BNavg,NBavg] = obj.cross_corr_helper(tetrode_cells,...
                        BBavg,NNavg,BNavg,NBavg);
                        tosave = 0;
                        if tosave
                            saveas(fig, xcorr_filename)
                            disp(xcorr_filename)
                        end
                        close(fig)

                    end
                end
            end
        end
        