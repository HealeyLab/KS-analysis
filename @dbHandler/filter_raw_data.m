function filter_raw_data(~)
           target_dir = pwd;
           
           [origFiles, origDataPath] = ... % crucial distinction: files vs file (current)
                uigetfile('*.rhd', 'Select an RHD2000 Data File', 'MultiSelect', 'on');
            cd(origDataPath)
            if iscell(origFiles)
                ts_label = origFiles{1};
            else
                ts_label = origFiles;
            end
            %%
            ts_label = ts_label(1:end-4);
            targetFileName = fullfile(target_dir, 'raw_filtered.dat'); 

            % open raw.dat to write
            fid = fopen(targetFileName, 'w'); % open .dat file for writing
            filearray = [];

            %ordering files
            for i = 1:length(origFiles)
                filearray = [filearray dir(char(origFiles{i}))];
            end

            [~, idx] = sort({filearray.date});
            origFiles = origFiles(idx);
            board_adc = [];
            for i=1:length(origFiles)
                [amplifier_data, frequency_parameters, board_adc_data] = read_Intan_RHD2000_file_MML_DJP(...
                    fullfile(filearray(i).folder, filearray(i).name),0);

                disp([num2str(i) ' of ' num2str(length(origFiles))])

                tic
                % only runs once
                if ~exist('a1', 'var')
                    [b1, a1] = butter(3, 300/frequency_parameters.amplifier_sample_rate*2, 'high');
                end
                % filter
                dataRAW = amplifier_data';
                dataRAW = single(dataRAW);

                datr = filter(b1, a1, dataRAW);
                datr = flipud(datr);
                datr = filter(b1, a1, datr);
                datr = flipud(datr);
                datr=datr';
                fwrite(fid, datr(:),'int16'); % append to .dat file
                %     fwrite(fid1a, amplifier_data(:),'int16'); % append to .dat file
                board_adc = [board_adc board_adc_data];
                toc
            end
            fclose(fid);
            cd(target_dir)
            beep
        end   