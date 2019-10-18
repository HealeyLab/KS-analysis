function gen_waveforms(obj, key)
keys = obj.get_key_family(key);

for i = 1:length(keys)    
    figure('units','normalized','outerposition',[0 0 .25 .8]);
    s = obj.db(keys{i});
    wf = s.spike_waveforms;
    wf_mean = nanmean(wf,2);
    hold on
   
    s.p2p = obj.get_p2p(wf_mean);
    s.sym = obj.get_sym(wf_mean);
    if s.p2p >= .43
        color = obj.BB;
    else
        color = obj.NN;
    end
    plot(wf, 'y')
    std_wf = nanstd(wf')';
    plot(wf_mean, 'LineWidth', 1.25, 'Color', color);
    plot(wf_mean+std_wf,'LineWidth', 1.25, 'Color', color);
    plot(wf_mean-std_wf,'LineWidth', 1.25, 'Color', color);  
    set(gca, 'XTick', [])
    ylim([(min(wf_mean-std_wf*2)) (max(wf_mean+std_wf)*2)]);
    
    key = replace(key, '_', ' ');
    title({key;['p2p: ' num2str(s.p2p) ';sym: ' num2str(s.sym)]});
end
end