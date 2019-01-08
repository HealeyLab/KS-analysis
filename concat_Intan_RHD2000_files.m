function [file, path] = concat_Intan_RHD2000_files

% concat_Intan_RHD2000_files
%
% I am blatantly editing this file, originally read_Intan_RHD2000_file.m,
% and making it my own. This file will concatenate multiple files by:
% 1) making multiselect 'on'
% 2) concatenating each instance of each and every variable that is saved
% in user's workspace. This will be done by creating a new workspace,
% concatenating them there, and porting the very final product to the
% workspace of the user.
% How? toward the end of the loop, a function (concat_struct) will be called that will
% concatenate each element to an element of a struct containing all
% attributes relevant.
% Once all the
% variables have been added, the constituent values (vectors or matrices
% will be concatenated and finally added to the 'base' workspace, or the
% workspace that the user can see.

[files, path] = ... % crucial distinction: files vs file (current)
    uigetfile('*.rhd', 'Select an RHD2000 Data File', 'MultiSelect', 'on');

if(iscell(files)) % when it's many files, it's number of files x 1
    loop_length = size(files,2);
elseif (ischar(files)) % when it's a single file, it's 1xlength of file string
    loop_length = 1;
else
    return;
end

notes_temp = [];
frequency_parameters_temp = [];
reference_channels_temp = [];
amplifier_channels_temp = [];
amplifier_data_temp = []; %
t_amplifier_temp = [];
spike_triggers_temp = struct(... % nested structure
    'voltage_trigger_mode', [],...
    'voltage_threshold', [],...
    'digital_trigger_channel', [],...
    'digital_edge_polarity', []);
aux_input_channels_temp = [];
aux_input_data_temp = [];
t_aux_input_temp = [];
supply_voltage_channels_temp = [];
supply_voltage_data_temp = [];
t_supply_voltage_temp = [];
board_adc_channels_temp = [];
board_adc_data_temp = [];
t_board_adc_temp = [];
board_dig_in_channels_temp = [];
board_dig_in_data_temp = [];
t_dig_temp = [];
board_dig_out_channels_temp = [];
board_dig_out_data_temp = [];
temp_sensor_data_temp = [];
t_temp_sensor_temp = [];

for index=1:loop_length
        
    % Read most recent file automatically.
    % path = 'C:\Users\Reid\Documents\RHD2132\testing\';
    % d = dir([path '*.rhd']);
    % file = d(end).name;
    tic;
    
    if (loop_length == 1) % if selected one, just an array
        file = files;
        filename = fullfile(path, file);
    else % if selected more than one, cells of arrays
        file = files(index);
        filename_cell = fullfile(path, file);
        filename = filename_cell{1}; % {} bc its a cell
    end
    
    fid = fopen(filename, 'r'); 
    
    s = dir(filename);
    filesize = s.bytes;
    
    % Check 'magic number' at beginning of file to make sure this is an Intan
    % Technologies RHD2000 data file.
    magic_number = fread(fid, 1, 'uint32');
    if magic_number ~= hex2dec('c6912702')
        error('Unrecognized file type.');
    end
    
    % Read version number.
    data_file_main_version_number = fread(fid, 1, 'int16');
    data_file_secondary_version_number = fread(fid, 1, 'int16');
    
    fprintf(1, '\n');
    fprintf(1, 'Reading Intan Technologies RHD2000 Data File, Version %d.%d\n', ...
        data_file_main_version_number, data_file_secondary_version_number);
    fprintf(1, '\n');
    
    if (data_file_main_version_number == 1)
        num_samples_per_data_block = 60;
    else
        num_samples_per_data_block = 128;
    end
    
    % Read information of sampling rate and amplifier frequency settings.
    sample_rate = fread(fid, 1, 'single');
    dsp_enabled = fread(fid, 1, 'int16');
    actual_dsp_cutoff_frequency = fread(fid, 1, 'single');
    actual_lower_bandwidth = fread(fid, 1, 'single');
    actual_upper_bandwidth = fread(fid, 1, 'single');
    
    desired_dsp_cutoff_frequency = fread(fid, 1, 'single');
    desired_lower_bandwidth = fread(fid, 1, 'single');
    desired_upper_bandwidth = fread(fid, 1, 'single');
    
    % This tells us if a software 50/60 Hz notch filter was enabled during
    % the data acquisition.
    notch_filter_mode = fread(fid, 1, 'int16');
    notch_filter_frequency = 0;
    if (notch_filter_mode == 1)
        notch_filter_frequency = 50;
    elseif (notch_filter_mode == 2)
        notch_filter_frequency = 60;
    end
    
    desired_impedance_test_frequency = fread(fid, 1, 'single');
    actual_impedance_test_frequency = fread(fid, 1, 'single');
    
    % Place notes in data strucure
    notes = struct( ...
        'note1', fread_QString(fid), ...
        'note2', fread_QString(fid), ...
        'note3', fread_QString(fid) );
    
    % If data file is from GUI v1.1 or later, see if temperature sensor data
    % was saved.
    num_temp_sensor_channels = 0;
    if ((data_file_main_version_number == 1 && data_file_secondary_version_number >= 1) ...
            || (data_file_main_version_number > 1))
        num_temp_sensor_channels = fread(fid, 1, 'int16');
    end
    
    % If data file is from GUI v1.3 or later, load eval board mode.
    eval_board_mode = 0;
    if ((data_file_main_version_number == 1 && data_file_secondary_version_number >= 3) ...
            || (data_file_main_version_number > 1))
        eval_board_mode = fread(fid, 1, 'int16');
    end
    
    % If data file is from v2.0 or later (Intan Recording Controller),
    % load name of digital reference channel.
    if (data_file_main_version_number > 1)
        reference_channel = fread_QString(fid);
    end
    
    % Place frequency-related information in data structure.
    frequency_parameters = struct( ...
        'amplifier_sample_rate', sample_rate, ...
        'aux_input_sample_rate', sample_rate / 4, ...
        'supply_voltage_sample_rate', sample_rate / num_samples_per_data_block, ...
        'board_adc_sample_rate', sample_rate, ...
        'board_dig_in_sample_rate', sample_rate, ...
        'desired_dsp_cutoff_frequency', desired_dsp_cutoff_frequency, ...
        'actual_dsp_cutoff_frequency', actual_dsp_cutoff_frequency, ...
        'dsp_enabled', dsp_enabled, ...
        'desired_lower_bandwidth', desired_lower_bandwidth, ...
        'actual_lower_bandwidth', actual_lower_bandwidth, ...
        'desired_upper_bandwidth', desired_upper_bandwidth, ...
        'actual_upper_bandwidth', actual_upper_bandwidth, ...
        'notch_filter_frequency', notch_filter_frequency, ...
        'desired_impedance_test_frequency', desired_impedance_test_frequency, ...
        'actual_impedance_test_frequency', actual_impedance_test_frequency );
    
    % Define data structure for spike trigger settings.
    spike_trigger_struct = struct( ...
        'voltage_trigger_mode', {}, ...
        'voltage_threshold', {}, ...
        'digital_trigger_channel', {}, ...
        'digital_edge_polarity', {} );
    
    new_trigger_channel = struct(spike_trigger_struct);
    spike_triggers = struct(spike_trigger_struct);
    
    % Define data structure for data channels.
    channel_struct = struct( ...
        'native_channel_name', {}, ...
        'custom_channel_name', {}, ...
        'native_order', {}, ...
        'custom_order', {}, ...
        'board_stream', {}, ...
        'chip_channel', {}, ...
        'port_name', {}, ...
        'port_prefix', {}, ...
        'port_number', {}, ...
        'electrode_impedance_magnitude', {}, ...
        'electrode_impedance_phase', {} );
    
    new_channel = struct(channel_struct);
    
    % Create structure arrays for each type of data channel.
    amplifier_channels = struct(channel_struct);
    aux_input_channels = struct(channel_struct);
    supply_voltage_channels = struct(channel_struct);
    board_adc_channels = struct(channel_struct);
    board_dig_in_channels = struct(channel_struct);
    board_dig_out_channels = struct(channel_struct);
    
    amplifier_index = 1;
    aux_input_index = 1;
    supply_voltage_index = 1;
    board_adc_index = 1;
    board_dig_in_index = 1;
    board_dig_out_index = 1;
    
    % Read signal summary from data file header.
    
    number_of_signal_groups = fread(fid, 1, 'int16');
    
    for signal_group = 1:number_of_signal_groups
        signal_group_name = fread_QString(fid);
        signal_group_prefix = fread_QString(fid);
        signal_group_enabled = fread(fid, 1, 'int16');
        signal_group_num_channels = fread(fid, 1, 'int16');
        signal_group_num_amp_channels = fread(fid, 1, 'int16');
        
        if (signal_group_num_channels > 0 && signal_group_enabled > 0)
            new_channel(1).port_name = signal_group_name;
            new_channel(1).port_prefix = signal_group_prefix;
            new_channel(1).port_number = signal_group;
            for signal_channel = 1:signal_group_num_channels
                new_channel(1).native_channel_name = fread_QString(fid);
                new_channel(1).custom_channel_name = fread_QString(fid);
                new_channel(1).native_order = fread(fid, 1, 'int16');
                new_channel(1).custom_order = fread(fid, 1, 'int16');
                signal_type = fread(fid, 1, 'int16');
                channel_enabled = fread(fid, 1, 'int16');
                new_channel(1).chip_channel = fread(fid, 1, 'int16');
                new_channel(1).board_stream = fread(fid, 1, 'int16');
                new_trigger_channel(1).voltage_trigger_mode = fread(fid, 1, 'int16');
                new_trigger_channel(1).voltage_threshold = fread(fid, 1, 'int16');
                new_trigger_channel(1).digital_trigger_channel = fread(fid, 1, 'int16');
                new_trigger_channel(1).digital_edge_polarity = fread(fid, 1, 'int16');
                new_channel(1).electrode_impedance_magnitude = fread(fid, 1, 'single');
                new_channel(1).electrode_impedance_phase = fread(fid, 1, 'single');
                
                if (channel_enabled)
                    switch (signal_type)
                        case 0
                            amplifier_channels(amplifier_index) = new_channel;
                            spike_triggers(amplifier_index) = new_trigger_channel;
                            amplifier_index = amplifier_index + 1;
                        case 1
                            aux_input_channels(aux_input_index) = new_channel;
                            aux_input_index = aux_input_index + 1;
                        case 2
                            supply_voltage_channels(supply_voltage_index) = new_channel;
                            supply_voltage_index = supply_voltage_index + 1;
                        case 3
                            board_adc_channels(board_adc_index) = new_channel;
                            board_adc_index = board_adc_index + 1;
                        case 4
                            board_dig_in_channels(board_dig_in_index) = new_channel;
                            board_dig_in_index = board_dig_in_index + 1;
                        case 5
                            board_dig_out_channels(board_dig_out_index) = new_channel;
                            board_dig_out_index = board_dig_out_index + 1;
                        otherwise
                            error('Unknown channel type');
                    end
                end
                
            end
        end
    end
    
    % Summarize contents of data file.
    num_amplifier_channels = amplifier_index - 1;
    num_aux_input_channels = aux_input_index - 1;
    num_supply_voltage_channels = supply_voltage_index - 1;
    num_board_adc_channels = board_adc_index - 1;
    num_board_dig_in_channels = board_dig_in_index - 1;
    num_board_dig_out_channels = board_dig_out_index - 1;
    
    fprintf(1, 'Found %d amplifier channel%s.\n', ...
        num_amplifier_channels, plural(num_amplifier_channels));
    fprintf(1, 'Found %d auxiliary input channel%s.\n', ...
        num_aux_input_channels, plural(num_aux_input_channels));
    fprintf(1, 'Found %d supply voltage channel%s.\n', ...
        num_supply_voltage_channels, plural(num_supply_voltage_channels));
    fprintf(1, 'Found %d board ADC channel%s.\n', ...
        num_board_adc_channels, plural(num_board_adc_channels));
    fprintf(1, 'Found %d board digital input channel%s.\n', ...
        num_board_dig_in_channels, plural(num_board_dig_in_channels));
    fprintf(1, 'Found %d board digital output channel%s.\n', ...
        num_board_dig_out_channels, plural(num_board_dig_out_channels));
    fprintf(1, 'Found %d temperature sensors channel%s.\n', ...
        num_temp_sensor_channels, plural(num_temp_sensor_channels));
    fprintf(1, '\n');
    
    % Determine how many samples the data file contains.
    
    % Each data block contains num_samples_per_data_block amplifier samples.
    bytes_per_block = num_samples_per_data_block * 4;  % timestamp data
    bytes_per_block = bytes_per_block + num_samples_per_data_block * 2 * num_amplifier_channels;
    % Auxiliary inputs are sampled 4x slower than amplifiers
    bytes_per_block = bytes_per_block + (num_samples_per_data_block / 4) * 2 * num_aux_input_channels;
    % Supply voltage is sampled once per data block
    bytes_per_block = bytes_per_block + 1 * 2 * num_supply_voltage_channels;
    % Board analog inputs are sampled at same rate as amplifiers
    bytes_per_block = bytes_per_block + num_samples_per_data_block * 2 * num_board_adc_channels;
    % Board digital inputs are sampled at same rate as amplifiers
    if (num_board_dig_in_channels > 0)
        bytes_per_block = bytes_per_block + num_samples_per_data_block * 2;
    end
    % Board digital outputs are sampled at same rate as amplifiers
    if (num_board_dig_out_channels > 0)
        bytes_per_block = bytes_per_block + num_samples_per_data_block * 2;
    end
    % Temp sensor is sampled once per data block
    if (num_temp_sensor_channels > 0)
        bytes_per_block = bytes_per_block + 1 * 2 * num_temp_sensor_channels;
    end
    
    % How many data blocks remain in this file?
    data_present = 0;
    bytes_remaining = filesize - ftell(fid);
    if (bytes_remaining > 0)
        data_present = 1;
    end
    
    num_data_blocks = bytes_remaining / bytes_per_block;
    
    num_amplifier_samples = num_samples_per_data_block * num_data_blocks;
    num_aux_input_samples = (num_samples_per_data_block / 4) * num_data_blocks;
    num_supply_voltage_samples = 1 * num_data_blocks;
    num_board_adc_samples = num_samples_per_data_block * num_data_blocks;
    num_board_dig_in_samples = num_samples_per_data_block * num_data_blocks;
    num_board_dig_out_samples = num_samples_per_data_block * num_data_blocks;
    
    record_time = num_amplifier_samples / sample_rate;
    
    if (data_present)
        fprintf(1, 'File contains %0.3f seconds of data.  Amplifiers were sampled at %0.2f kS/s.\n', ...
            record_time, sample_rate / 1000);
        fprintf(1, '\n');
    else
        fprintf(1, 'Header file contains no data.  Amplifiers were sampled at %0.2f kS/s.\n', ...
            sample_rate / 1000);
        fprintf(1, '\n');
    end
    
    if (data_present)
        
        % Pre-allocate memory for data.
        fprintf(1, 'Allocating memory for data...\n');
        
        t_amplifier = zeros(1, num_amplifier_samples);
        
        amplifier_data = zeros(num_amplifier_channels, num_amplifier_samples);
        aux_input_data = zeros(num_aux_input_channels, num_aux_input_samples);
        supply_voltage_data = zeros(num_supply_voltage_channels, num_supply_voltage_samples);
        temp_sensor_data = zeros(num_temp_sensor_channels, num_supply_voltage_samples);
        board_adc_data = zeros(num_board_adc_channels, num_board_adc_samples);
        board_dig_in_data = zeros(num_board_dig_in_channels, num_board_dig_in_samples);
        board_dig_in_raw = zeros(1, num_board_dig_in_samples);
        board_dig_out_data = zeros(num_board_dig_out_channels, num_board_dig_out_samples);
        board_dig_out_raw = zeros(1, num_board_dig_out_samples);
        
        % Read sampled data from file.
        fprintf(1, 'Reading data from file...\n');
        
        amplifier_index = 1;
        aux_input_index = 1;
        supply_voltage_index = 1;
        board_adc_index = 1;
        board_dig_in_index = 1;
        board_dig_out_index = 1;
        
        print_increment = 10;
        percent_done = print_increment;
        for i=1:num_data_blocks
            % In version 1.2, we moved from saving timestamps as unsigned
            % integeters to signed integers to accomidate negative (adjusted)
            % timestamps for pretrigger data.
            if ((data_file_main_version_number == 1 && data_file_secondary_version_number >= 2) ...
                    || (data_file_main_version_number > 1))
                t_amplifier(amplifier_index:(amplifier_index + num_samples_per_data_block - 1)) = fread(fid, num_samples_per_data_block, 'int32');
            else
                t_amplifier(amplifier_index:(amplifier_index + num_samples_per_data_block - 1)) = fread(fid, num_samples_per_data_block, 'uint32');
            end
            if (num_amplifier_channels > 0)
                amplifier_data(:, amplifier_index:(amplifier_index + num_samples_per_data_block - 1)) = fread(fid, [num_samples_per_data_block, num_amplifier_channels], 'uint16')';
            end
            if (num_aux_input_channels > 0)
                aux_input_data(:, aux_input_index:(aux_input_index + (num_samples_per_data_block / 4) - 1)) = fread(fid, [(num_samples_per_data_block / 4), num_aux_input_channels], 'uint16')';
            end
            if (num_supply_voltage_channels > 0)
                supply_voltage_data(:, supply_voltage_index) = fread(fid, [1, num_supply_voltage_channels], 'uint16')';
            end
            if (num_temp_sensor_channels > 0)
                temp_sensor_data(:, supply_voltage_index) = fread(fid, [1, num_temp_sensor_channels], 'int16')';
            end
            if (num_board_adc_channels > 0)
                board_adc_data(:, board_adc_index:(board_adc_index + num_samples_per_data_block - 1)) = fread(fid, [num_samples_per_data_block, num_board_adc_channels], 'uint16')';
            end
            if (num_board_dig_in_channels > 0)
                board_dig_in_raw(board_dig_in_index:(board_dig_in_index + num_samples_per_data_block - 1)) = fread(fid, num_samples_per_data_block, 'uint16');
            end
            if (num_board_dig_out_channels > 0)
                board_dig_out_raw(board_dig_out_index:(board_dig_out_index + num_samples_per_data_block - 1)) = fread(fid, num_samples_per_data_block, 'uint16');
            end
            
            amplifier_index = amplifier_index + num_samples_per_data_block;
            aux_input_index = aux_input_index + (num_samples_per_data_block / 4);
            supply_voltage_index = supply_voltage_index + 1;
            board_adc_index = board_adc_index + num_samples_per_data_block;
            board_dig_in_index = board_dig_in_index + num_samples_per_data_block;
            board_dig_out_index = board_dig_out_index + num_samples_per_data_block;
            
            fraction_done = 100 * (i / num_data_blocks);
            if (fraction_done >= percent_done)
                fprintf(1, '%d%% done...\n', percent_done);
                percent_done = percent_done + print_increment;
            end
        end
        
        % Make sure we have read exactly the right amount of data.
        bytes_remaining = filesize - ftell(fid);
        if (bytes_remaining ~= 0)
            %error('Error: End of file not reached.');
        end
        
    end
    
    % Close data file.
    fclose(fid);
    
    if (data_present)
        
        fprintf(1, 'Parsing data...\n');
        
        % Extract digital input channels to separate variables.
        for i=1:num_board_dig_in_channels
            mask = 2^(board_dig_in_channels(i).native_order) * ones(size(board_dig_in_raw));
            board_dig_in_data(i, :) = (bitand(board_dig_in_raw, mask) > 0);
        end
        for i=1:num_board_dig_out_channels
            mask = 2^(board_dig_out_channels(i).native_order) * ones(size(board_dig_out_raw));
            board_dig_out_data(i, :) = (bitand(board_dig_out_raw, mask) > 0);
        end
        
        % Scale voltage levels appropriately.
        amplifier_data = 0.195 * (amplifier_data - 32768); % units = microvolts
        aux_input_data = 37.4e-6 * aux_input_data; % units = volts
        supply_voltage_data = 74.8e-6 * supply_voltage_data; % units = volts
        if (eval_board_mode == 1)
            board_adc_data = 152.59e-6 * (board_adc_data - 32768); % units = volts
        elseif (eval_board_mode == 13) % Intan Recording Controller
            board_adc_data = 312.5e-6 * (board_adc_data - 32768); % units = volts
        else
            board_adc_data = 50.354e-6 * board_adc_data; % units = volts
        end
        temp_sensor_data = temp_sensor_data / 100; % units = deg C
        
        % Check for gaps in timestamps.
        num_gaps = sum(diff(t_amplifier) ~= 1);
        if (num_gaps == 0)
            fprintf(1, 'No missing timestamps in data.\n');
        else
            fprintf(1, 'Warning: %d gaps in timestamp data found.  Time scale will not be uniform!\n', ...
                num_gaps);
        end
        
        % Scale time steps (units = seconds).
        t_amplifier = t_amplifier / sample_rate;
        t_aux_input = t_amplifier(1:4:end);
        t_supply_voltage = t_amplifier(1:num_samples_per_data_block:end);
        t_board_adc = t_amplifier;
        t_dig = t_amplifier;
        t_temp_sensor = t_supply_voltage;
        
        % If the software notch filter was selected during the recording, apply the
        % same notch filter to amplifier data here.
        if (notch_filter_frequency > 0)
            fprintf(1, 'Applying notch filter...\n');
            
            print_increment = 10;
            percent_done = print_increment;
            for i=1:num_amplifier_channels
                amplifier_data(i,:) = ...
                    notch_filter(amplifier_data(i,:), sample_rate, notch_filter_frequency, 10);
                
                fraction_done = 100 * (i / num_amplifier_channels);
                if (fraction_done >= percent_done)
                    fprintf(1, '%d%% done...\n', percent_done);
                    percent_done = percent_done + print_increment;
                end
                
            end
        end
        
    end
    %% Concat variables in struct.
    if exist('notes', 'var')
        notes_temp = [notes_temp; notes]; % creates length(files) x 1 struct
        disp 'notes, '
    end
    if exist('frequency_parameters', 'var')
        frequency_parameters_temp = frequency_parameters; % no concat
        disp 'frequency params, '
    end
    if exist('reference_channel', 'var')
        reference_channel_temp = reference_channel; % no concat
        disp 'ref chan, '
    end
    if exist('amplifier_channels', 'var')
        amplifier_channels_temp = amplifier_channels; % no concat
        disp 'amp chan, '
    end
    if exist('amplifier_data', 'var')
        amplifier_data_temp = [amplifier_data_temp amplifier_data]; % actually horizontal
        disp 'amp data, '
    end
    if exist('t_amplifier', 'var')% also horizontal, is just time indices
        t_amplifier_temp = time_indexes(t_amplifier_temp, t_amplifier, index);     
        disp 't_amp, '
    end
    if exist('spike_triggers', 'var') % this seems to actually concatenate the values of each field
        spike_triggers_temp.voltage_trigger_mode =...
            [spike_triggers_temp.voltage_trigger_mode spike_triggers.voltage_trigger_mode];
        
        spike_triggers_temp.voltage_threshold =...
            [spike_triggers_temp.voltage_threshold spike_triggers.voltage_threshold];
        
        spike_triggers_temp.digital_trigger_channel =...
            [spike_triggers_temp.digital_trigger_channel spike_triggers.digital_trigger_channel];
        
        spike_triggers_temp.digital_edge_polarity =...
            [spike_triggers_temp.digital_edge_polarity spike_triggers.digital_edge_polarity];
        disp 'spike trig, '
    end
    if exist('aux_input_channels', 'var')
        aux_input_channels_temp = aux_input_channels;
        disp 'aux in chans, '
    end
    if exist('aux_input_data', 'var')
        aux_input_data_temp = [aux_input_data_temp aux_input_data];
        disp 'aux in data, '
    end
    if exist('t_aux_input', 'var')
        t_aux_input_temp = time_indexes(t_aux_input_temp, t_aux_input, index); %variables_to_concat.t_aux_input = [variables_to_concat.t_aux_input t_aux_input];
        disp 't aux input, '
    end
    if exist('supply_voltage_channels', 'var')
        supply_voltage_channels_temp = supply_voltage_channels;
        disp 'supply volt chans, '
    end 
    if exist('supply_voltage_data', 'var')
        supply_voltage_data_temp = [supply_voltage_data_temp supply_voltage_data];
        disp 'supply volt chans, '
    end 
    if exist('t_supply_voltage', 'var')
        t_supply_voltage_temp = time_indexes(t_supply_voltage_temp, t_supply_voltage, index);
        disp 't supply volt, '
    end 
    if exist('board_adc_channels', 'var')
        board_adc_channels_temp = board_adc_channels;
        disp 'board adc chans, '
    end 
    if exist('board_adc_data', 'var') % this is where the money happens
        board_adc_data_temp =...
            [board_adc_data_temp board_adc_data];
        disp ' board adc data, '
    end 
    if exist('t_board_adc', 'var')
        t_board_adc_temp = time_indexes(t_board_adc_temp, t_board_adc, index);
        disp ' t board adc, '
    end 
    if exist('board_dig_in_channels', 'var')
        board_dig_in_channels_temp =...
            [board_dig_in_channels_temp board_dig_in_channels];
        disp 'board dig in chans, '
    end 
    if exist('board_dig_in_data', 'var')
        board_dig_in_data_temp =...
            [board_dig_in_data_temp board_dig_in_data];
        disp 'board dig in data, '
    end 
    if exist('t_dig', 'var')
        t_dig_temp = time_indexes(t_dig_temp, t_dig, index);       
        disp 't dig, '
    end
    if exist('board_dig_out_channels', 'var')
        board_dig_in_channels_temp =...
            [board_dig_in_channels_temp board_dig_in_channels];
        disp 'board dig out chans, '
    end 
    if exist('board_dig_out_data', 'var')
        board_dig_in_data_temp =...
            [board_dig_in_data_temp board_dig_in_data];
        disp 'board dig out data, '
    end 
    if exist('temp_sensor_data', 'var')
        temp_sensor_data_temp =...
            [temp_sensor_data_temp temp_sensor_data];
        disp 'temp sensor data, '
    end
    if exist('t_temp_sensor', 'var')
        t_temp_sensor_temp = time_indexes(t_temp_sensor_temp, t_temp_sensor, index);
        disp 't temp sensor, '
    end
end % end of for loop through all files selected


%% Move variables to base workspace.
% Changed variable names to var_temp
notes = notes_temp;
frequency_parameters = frequency_parameters_temp;
reference_channels = reference_channels_temp;
amplifier_channels = amplifier_channels_temp;
amplifier_data = amplifier_data_temp;
t_amplifier = t_amplifier_temp;
spike_triggers = spike_triggers_temp;
aux_input_channels = aux_input_channels_temp;
aux_input_data = aux_input_data_temp;
t_aux_input = t_aux_input_temp;
supply_voltage_channels = supply_voltage_channels_temp;
supply_voltage_data = supply_voltage_data_temp;
t_supply_voltage = t_supply_voltage_temp; 
board_adc_channels = board_adc_channels_temp;
board_adc_data = board_adc_data_temp;
t_board_adc = t_board_adc_temp;
board_dig_in_channels = board_dig_in_channels_temp;
board_dig_in_data = board_dig_in_data_temp;
t_dig = t_dig_temp;
board_dig_out_channels = board_dig_out_channels_temp;
board_dig_out_data = board_dig_out_data_temp;
temp_sensor_data = temp_sensor_data_temp;
t_temp_sensor = t_temp_sensor_temp;

move_to_base_workspace(notes);
move_to_base_workspace(frequency_parameters);
if (data_file_main_version_number > 1)
    move_to_base_workspace(reference_channel);
end

if (num_amplifier_channels > 0)
    move_to_base_workspace(amplifier_channels);
    if (data_present)
        move_to_base_workspace(amplifier_data);
        move_to_base_workspace(t_amplifier);
    end
    move_to_base_workspace(spike_triggers);
end
if (num_aux_input_channels > 0)
    move_to_base_workspace(aux_input_channels);
    if (data_present)
        move_to_base_workspace(aux_input_data);
        move_to_base_workspace(t_aux_input);
    end
end
if (num_supply_voltage_channels > 0)
    move_to_base_workspace(supply_voltage_channels);
    if (data_present)
        move_to_base_workspace(supply_voltage_data);
        move_to_base_workspace(t_supply_voltage);
    end
end
if (num_board_adc_channels > 0)
    move_to_base_workspace(board_adc_channels);
    if (data_present)
        move_to_base_workspace(board_adc_data);
        move_to_base_workspace(t_board_adc);
    end
end
if (num_board_dig_in_channels > 0)
    move_to_base_workspace(board_dig_in_channels);
    if (data_present)
        move_to_base_workspace(board_dig_in_data);
        move_to_base_workspace(t_dig);
    end
end
if (num_board_dig_out_channels > 0)
    move_to_base_workspace(board_dig_out_channels);
    if (data_present)
        move_to_base_workspace(board_dig_out_data);
        move_to_base_workspace(t_dig);
    end
end
if (num_temp_sensor_channels > 0)
    if (data_present)
        move_to_base_workspace(temp_sensor_data);
        move_to_base_workspace(t_temp_sensor);
    end
end

clear *_temp % get rid of temp variables

fprintf(1, 'Done!  Elapsed time: %0.1f seconds\n', toc);
if (data_present)
    fprintf(1, 'Extracted data are now available in the MATLAB workspace.\n');
else
    fprintf(1, 'Extracted waveform information is now available in the MATLAB workspace.\n');
end
fprintf(1, 'Type ''whos'' to see variables.\n');
fprintf(1, '\n');

return


function a = fread_QString(fid)

% a = read_QString(fid)
%
% Read Qt style QString.  The first 32-bit unsigned number indicates
% the length of the string (in bytes).  If this number equals 0xFFFFFFFF,
% the string is null.

a = '';
length = fread(fid, 1, 'uint32');
if length == hex2num('ffffffff')
    return;
end
% convert length from bytes to 16-bit Unicode words
length = length / 2;

for i=1:length
    a(i) = fread(fid, 1, 'uint16');
end

return


function s = plural(n)

% s = plural(n)
%
% Utility function to optionally plurailze words based on the value
% of n.

if (n == 1)
    s = '';
else
    s = 's';
end

return


function out = notch_filter(in, fSample, fNotch, Bandwidth)

% out = notch_filter(in, fSample, fNotch, Bandwidth)
%
% Implements a notch filter (e.g., for 50 or 60 Hz) on vector 'in'.
% fSample = sample rate of data (in Hz or Samples/sec)
% fNotch = filter notch frequency (in Hz)
% Bandwidth = notch 3-dB bandwidth (in Hz).  A bandwidth of 10 Hz is
%   recommended for 50 or 60 Hz notch filters; narrower bandwidths lead to
%   poor time-domain properties with an extended ringing response to
%   transient disturbances.
%
% Example:  If neural data was sampled at 30 kSamples/sec
% and you wish to implement a 60 Hz notch filter:
%
% out = notch_filter(in, 30000, 60, 10);

tstep = 1/fSample;
Fc = fNotch*tstep;

L = length(in);

% Calculate IIR filter parameters
d = exp(-2*pi*(Bandwidth/2)*tstep);
b = (1 + d*d)*cos(2*pi*Fc);
a0 = 1;
a1 = -b;
a2 = d*d;
a = (1 + d*d)/2;
b0 = 1;
b1 = -2*cos(2*pi*Fc);
b2 = 1;

out = zeros(size(in));
out(1) = in(1);
out(2) = in(2);
% (If filtering a continuous data stream, change out(1) and out(2) to the
%  previous final two values of out.)

% Run filter
for i=3:L
    out(i) = (a*b2*in(i-2) + a*b1*in(i-1) + a*b0*in(i) - a2*out(i-2) - a1*out(i-1))/a0;
end

return


function variable_name = get_variable_name(variable)
variable_name = inputname(1);
return;


function move_to_base_workspace(variable)

% move_to_base_workspace(variable)
%
% Move variable from function workspace to base MATLAB workspace so
% user will have access to it after the program ends.

variable_name = inputname(1);
assignin('base', variable_name, variable);
return;


function to_concat = time_indexes(to_concat, t_var, index)
    if (index == 1) % also horizontal, is just time indices
        to_concat = [to_concat t_var]; 
    elseif(index == 2) % also horizontal, is just time indices
        to_concat = [to_concat to_concat(end)+t_var]; 
    end
return