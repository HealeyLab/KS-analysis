function gen_waveforms(obj, key)
keys = obj.get_key_family(key);
for i = l:length(keys)    
    figure;
    s = obj.db(s);
    wf = s.spike_waveforms;
    wf_mean = nanmean(wf,2);
    hold on
    if s.p2p >= .43
        color = BB;
    else
        color = NN;
    end
    plot(wf_mean, 'LineWidth', 1.25, 'Color', color);
    plot(wf_mean+nanstd(wf')','LineWidth', 1.25, 'Color', color);
    plot(wf_mean-nanstd(wf')','LineWidth', 1.25, 'Color', color);  
    set(gca, 'XTick', [])
    title(['p2p: ' num2str(s.p2p) '            sym: ' num2str(s.sym)])
end
end