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
    key = keys{i};
    obj.plot_waveform(key);

    s = obj.db(key);
    key_str = replace(key, '_', ' ');
    title({key_str;['p2p: ' num2str(obj.get_p2p(s)) ';sym: ' num2str(obj.get_sym(s))]});
end
end