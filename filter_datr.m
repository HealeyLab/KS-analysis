function datr_out = filter_datr(amplifier_data, frequency_parameters)
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
    datr_out=datr';
    
end