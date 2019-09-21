function output = get_audio(~, aPath)
    wav_files = dir([aPath '\*.wav']);
    numfiles = length(wav_files);
    output = struct('wav', cell(numfiles,1),...
        'fs', cell(numfiles,1),...
        'name', cell(numfiles,1));
    for i = 1:numfiles
        [cur_wav, fs] = audioread(fullfile(wav_files(i).folder,...
            wav_files(i).name));

        output(i).wav = cur_wav(:,1);
        output(i).fs = fs;
        output(i).name = wav_files(i).name;
    end
end