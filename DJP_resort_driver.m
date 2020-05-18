%% Author: DJP
% I could have sworn I wrote this script already, but here I am, writing it
% again. Jeez. Gosh.
% This script is intended for re-running sorting on files and setting the
% settings to the chic settings du jour
pwd_split = split(pwd, '_');    
%     {'E:\DJP thesis sorting\mdk\3 4 20\mdk  3 4 30 habitu'}
%     {'200304'                                             }
%     {'164355'                                             }
%     {'Kilosort'                                           }

% get the timestamp key from filename
ts_key = pwd_split{end - 1}; % '164355'

file = dir(['../*' ts_key '*.rhd']);

[amplifier_data, ~, ~ ] = read_Intan_RHD2000_file_MML_DJP(fullfile(file.folder, file.name),0);
% std_fig = visualize_std(filter_datr(amplifier_data, frequency_parameters));

%% change threshold in COnfig file
fid = fopen('StandardConfig_MOVEME.m','r');
i = 1;
tline = fgetl(fid);
A{i} = tline;
while ischar(tline)
    i = i+1;
    tline = fgetl(fid);
    A{i} = tline;
end
fclose(fid);

A{7} = sprintf('%s',['ops.fbinary = ''' fullfile(pwd, 'raw_filtered.dat'';')]);
A{8} = sprintf('%s',['ops.fproc = ''' fullfile(pwd, 'temp_wh.dat'';')]);
A{9} = sprintf('%s', ['ops.root = ''' pwd ''';']);
A{24} = sprintf('%s',['ops.chanMap = ''' fullfile(pwd, 'chanMap.mat') ''';']);
% because we are doing initialize = fromData, this is to increase merging
A{14} = sprintf('ops.Nfilt = 32;');

% disp('close the fig')
% waitfor(std_fig)
% thres = input('what do you want the threshold to be?\n','s');
% if strcmp(thres,"")
%     thres = -6;
% end
A{47} = sprintf('ops.mergeT           = .08;');
A{48} = sprintf('ops.splitT           = .12;');
A{51} = sprintf('ops.initialize      = ''fromData'';');
A{52} = sprintf('ops.spkTh           = -5;');
% A{52} = sprintf('ops.spkTh           = -%s;      ', thres);

% Write cell A into txt
fid = fopen('StandardConfig_MOVEME.m', 'w');
for i = 1:numel(A)
    if A{i+1} == -1
        fprintf(fid,'%s', A{i});
        break
    else
        fprintf(fid,'%s\n', A{i});
    end
end
fclose(fid);
%% Edit masterfile
%% Third, edicd t the master file
fid5 = fopen('master_file_example_MOVEME.m','r');
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

C{6} = sprintf('pathToYourConfigFile = ''%s''; ', pwd);

% Write cell C into txt
fid5 = fopen('master_file_example_MOVEME.m', 'w');
for i = 1:numel(C)
    if C{i+1} == -1
        fprintf(fid5,'%s', C{i});
        break
    else
        fprintf(fid5,'%s\n', C{i});
    end
end
fclose(fid5);
%% close all
fclose('all');

%% run
run('master_file_example_MOVEME.m');
pushBulletDriver(strjoin(['done sorting ' string(pwd)]));