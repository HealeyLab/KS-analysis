function figs_out = gen_waveforms(obj, key, only)

if only
    keys = {key};
else
    keys = obj.get_key_family(key);
end

figs_out = [];
for i = 1:length(keys)    
    fig = figure('units','normalized','outerposition',[0 0 .25 .8]);
    figs_out = [figs_out fig];
    s = obj.db(keys{i});
    wf = s.spike_waveforms;
    wf_mean = nanmean(wf,2);
    hold on
   
    s.p2p = obj.get_p2p(s);
    s.sym = obj.get_sym(s);
    
    plot(wf, 'y', 'Color', [0,0,0]+.75)
    wf_std = nanstd(wf')';
    plot(wf_mean, 'LineWidth', 1.25, 'Color', 'k');
    plot(wf_mean+wf_std,'LineWidth', 1.25, 'Color', 'k');
    plot(wf_mean-wf_std,'LineWidth', 1.25, 'Color', 'k');  
    set(gca, 'XTick', [])
    ylim([(min(wf_mean-wf_std*5)) (max(wf_mean+wf_std)*5)]);
    
    key = replace(key, '_', ' ');
    title({key;['p2p: ' num2str(s.p2p) ';sym: ' num2str(s.sym)]});
    
end
end