function showSpectrogram(~, mic, adc_sr)
    
%    https://stackoverflow.com/questions/56416877/changing-axis-units-of-matlab-spectrogram
    [~,F,T,P] = spectrogram( mic, 256, [],[], adc_sr, 'yaxis');
    imagesc(T, F, 10*log10(P+eps)); % add eps like pspectrogram does
    axis xy
    ylabel('Frequency (Hz)')
    xlabel('Time (s)')
    h = colorbar;
    h.Label.String = 'Power/frequency (dB/Hz)';
    colorbar('delete');
end