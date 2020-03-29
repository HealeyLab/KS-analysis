function song_filtered = filter_song(~, board_adc, sr)
%% FILTER_SONG Filters song trace 300-12k Hz
% For visualizing purposes, I often need to clean up a song trace.

    mic = board_adc(end,:); % if you re trying to access the audio trace, do (1,:)
    fcutlow = 300;
    fcuthigh = 12e3;
    order = 3;
    try
        [b,a]=butter(order,[fcutlow,fcuthigh]/(sr/2),'bandpass');
    catch % sampling rate is too low for 12e3 cutoff
        [b,a]=butter(order,[fcutlow,9.99e3]/(sr/2),'bandpass');
    end
    mic = filter(b, a, mic);
    mic = fliplr(mic);
    mic = filter(b, a, mic);
    song_filtered = fliplr(mic);

end