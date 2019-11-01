function [fig_out] = gen_strf(obj, key_pattern)

keys = obj.show_keys(key_pattern);
for key_inds = 1:length(keys)
    key = keys{key_inds};
    s = obj.db(key);

    [on, off, wav_files] = obj.get_stim_on_off(key);

    % for each class of stim:
    si = s.stim_identities{1};
    uniq = unique(si); % remove redundancies
    usi = uniq(~contains(uniq, 'wav')); % remove non-non-wav files
    Z = nan(7,39); 
    % for each amp/freq permutation:
    for i = 1:length(usi)
        % for each stim itself
        ind = find(strcmp(si,usi{i}));

        sp_ts = s.spike_timestamps;
        start = on(ind);
        stop = off(ind);

        raster_arr = ((intersect(...
            sp_ts(sp_ts > start),...
            sp_ts(sp_ts < stop))...
            - start) / s.amplifier_sampling_rate)';                        

        % get x and y indices
        s_text = usi{i};
        spl = split(s_text, ' ');
        freq = str2double(spl{1});
        amp = str2double(spl{2});
        col = int16(log2(freq/200)/.1);
        row = 8-int16(amp / 10);

        Z(row,col)= length(raster_arr);

    end

    fig_out = figure('units','normalized','outerposition',[0 0 .25 .8]);

    xvar = int16(200 .* 2 .^ (.1* [1:39])); yvar = fliplr(int16(single(10 * [1:7])));

    color = obj.get_color(s);
    
    cc = max(max(Z));
    heatmap(xvar,yvar,Z, 'Colormap', fliplr([0:1/cc:1; 0:1/cc:1;0:1/cc:1] .* color')');

    xlabel('Frequency (Hz)')
    key = replace(key, '_', ' ');
    title({key});
end
end

