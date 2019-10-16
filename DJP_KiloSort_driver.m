%% Author: DJP
close all;

%% Conversion
[files, origDataPath] = ... % crucial distinction: files vs file (current)
    uigetfile('*.rhd', 'Select an RHD2000 Data File', 'MultiSelect', 'on');
% [file,origDataPath] = concat_Intan_RHD2000_files;
% make folder containing that files
if iscell(files)
    ts_label = files{1};
else
    ts_label = files;
end
ts_label = ts_label(1:end-4);
dataPath = fullfile(origDataPath,[ts_label '_Kilosort']);
mkdir(dataPath)
addpath(dataPath)
%%
dataFileName = fullfile(dataPath, 'raw_filtered.dat'); % [dataPath '\Kilosort_alldata\raw.dat'] % make .dat file
% saving adc data, too
board_adc = [];
% open raw.dat to write
fid = fopen(dataFileName, 'w'); % open .dat file for writing

filearray = [];
for i = 1:length(files)
    filearray = [filearray dir(char(files(i)))];
end

[~, idx] = sort({filearray.date});
files = files(idx);
for i=1:length(files)
    read_Intan_RHD2000_file_MML_DJP(fullfile(filearray(i).folder,files{i}),0);
%     
%     % only runs once
%     if ~exist('a1', 'var')
%         [b1, a1] = butter(3, 300/frequency_parameters.amplifier_sample_rate*2, 'high');
%     end
% 
%     % filter
%     dataRAW = amplifier_data';
%     dataRAW = single(dataRAW);
% 
%     datr = filter(b1, a1, dataRAW);
%     datr = flipud(datr);
%     datr = filter(b1, a1, datr);
%     datr = flipud(datr);
%     datr=datr';
%     fwrite(fid, datr(:),'int16'); % append to .dat file
    %     fwrite(fid1a, amplifier_data(:),'int16'); % append to .dat file
    board_adc = [board_adc board_adc_data];
end

fclose(fid);

adc_sr=frequency_parameters.board_adc_sample_rate;
save(fullfile(dataPath, 'adc_data'), 'board_adc', 'adc_sr', '-v7.3')
%% make params.py
fid2 = fopen(fullfile(dataPath, 'params.py'),'w');

dat_path = 'raw_filtered.dat';
n_channels_dat = string(length(amplifier_channels));
dtype = 'int16';
offset = 0;
sample_rate = string(frequency_parameters.amplifier_sample_rate);
hp_filtered = 'True';

fprintf(fid2, 'dat_path = r''%s''\n', dat_path);
fprintf(fid2, 'n_channels_dat = %s\n', n_channels_dat); %2d means two digit
fprintf(fid2, 'dtype = ''%s''\n', dtype);
fprintf(fid2, 'offset = 0\n');
fprintf(fid2, 'sample_rate = %s.\n', sample_rate);
fprintf(fid2, 'hp_filtered = %s\n', hp_filtered);

fclose(fid2);
clear fid2 dat_path dtype offset hp_filtered 
clear amplifier_channels amplifier_data aux_input_channels aux_input_data ...
        board_dig_in_data board_dig_in_channels filename frequency_parameters ...
        notes reference_channel spike_triggers supply_voltage_channels supply_voltage_data ...
        t_amplifier t_aux_input t_dig t_supply_voltage
%% Run Kilosort
% copy master file example and  standard config and then edit them
working_dir = 'C:\Users\danpo\Documents\MATLAB\DJP_KiloSort\KS-analysis';

ChannelMapFile_orig      = fullfile(working_dir, 'createChannelMapFile.m');
master_file_example_orig = fullfile(working_dir, 'master_file_example_MOVEME.m');
StandardConfig_orig      = fullfile(working_dir, 'StandardConfig_MOVEME.m');

copyfile(ChannelMapFile_orig,       dataPath)
copyfile(master_file_example_orig,  dataPath)
copyfile(StandardConfig_orig,       dataPath)

ChannelMapFile_pasted      = fullfile(dataPath, 'createChannelMapFile.m');
master_file_example_pasted = fullfile(dataPath, 'master_file_example_MOVEME.m');
StandardConfig_pasted      = fullfile(dataPath, 'StandardConfig_MOVEME.m');

%% First, StandardConfig is copied into cell with code taken from online
% https://www.mathworks.com/matlabcentral/answers/62986-how-to-change-a-specific-line-in-a-text-file
fid3 = fopen(StandardConfig_pasted,'r');
i = 1;
tline = fgetl(fid3);
A{i} = tline;
while ischar(tline)
    i = i+1;
    tline = fgetl(fid3);
    A{i} = tline;
end
fclose(fid3);
% Change lines of .m file
A{1}  = sprintf('ops.GPU=0');
A{7}  = sprintf('ops.fbinary = ''%s'';', dataFileName); % dataFileName = 'C:\Users\danpo\Documents\sorting\mdx 10 8 18\filename'
A{8}  = sprintf('ops.fproc = ''%s'';', fullfile(dataPath,'temp_wh.dat')); % dataPath = 'C:\Users\danpo\Documents\sorting\mdx 10 8 18\'
A{9}  = sprintf('ops.root = ''%s'';',dataPath);
A{11} = sprintf('ops.fs = %s;',sample_rate);
% change 12-14 to nothing if already covered in other thing
% NOTE: ONLY 16 CHANNELS FOR THESE RECORDINGS
A{12} = sprintf('ops.NchanTOT = %s;', n_channels_dat);
A{13} = sprintf('ops.Nchan = %s;', string(16));
A{14} = sprintf('ops.Nfilt = %s;',  string(str2num(n_channels_dat) * 2));% number of clusters to use (2-4 times more than Nchan, should be a multiple of 32)     		

A{24} = sprintf('ops.chanMap = ''%s'';',  fullfile(dataPath, 'chanMap.mat'));

% Write cell A into txt
fid3 = fopen(StandardConfig_pasted, 'w');
for i = 1:numel(A)
    if A{i+1} == -1
        fprintf(fid3,'%s', A{i});
        break
    else
        fprintf(fid3,'%s\n', A{i});
    end
end
fclose(fid3);
%% Second, ChannelMap is copied into cell with code taken from online
fid4 = fopen(ChannelMapFile_pasted,'r');
i = 1;
tline = fgetl(fid4);
B{i} = tline;
while ischar(tline)
    i = i+1;
    tline = fgetl(fid4);
    B{i} = tline;
end
fclose(fid4);
% Change lines of .m file
B{6} = sprintf('pathToYourConfigFile = ''%s'';', dataPath);
B{11} = sprintf('fs = %s;', sample_rate);
B{48} = sprintf('save(''%s'', ...',  fullfile(dataPath, 'chanMap.mat'));

% Write cell B into .txt
fid4 = fopen(ChannelMapFile_pasted, 'w');
for i = 1:numel(B)
    if B{i+1} == -1
        fprintf(fid4,'%s', B{i});
        break
    else
        fprintf(fid4,'%s\n', B{i});
    end
end
fclose(fid4);
%% Third, edit the master file
fid5 = fopen(master_file_example_pasted,'r');
i = 1;
tline = fgetl(fid5);
C{i} = tline;
while ischar(tline)
    i = i+1;
    tline = fgetl(fid5);
    C{i} = tline;
end
fclose(fid5);
% Change lines of .m file
C{6} = sprintf('pathToYourConfigFile = ''%s''; ', dataPath);
% C{24} = sprintf('rez = merge_posthoc2(rez);', dataPath);

% Write cell C into txt
fid5 = fopen(master_file_example_pasted, 'w');
for i = 1:numel(C)
    if C{i+1} == -1
        fprintf(fid5,'%s', C{i});
        break
    else
        fprintf(fid5,'%s\n', C{i});
    end
end
fclose(fid5);
%%
fclose('all');
%% Go
run(ChannelMapFile_pasted)
master_file_example_MOVEME
beep
pushBulletDriver('done');