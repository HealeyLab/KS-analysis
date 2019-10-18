function gen_strf(obj, key)

s = obj.db(key);

[on, off, wav_files] = obj.get_stim_on_off(key);

% for each class of stim:
si = s.stim_identities{1};
uniq = unique(si); % remove redundancies
usi = uniq(~contains(uniq, 'wav')); % remove non-non-wav files
Z = nan(39,10); % 39 columns (x) by 10 rows (y)
% for each amp/freq permutation:
for i = 1:length(usi)
    % for each stim itself
    ind = find(strcmp(si,usi{i}));
    
    % for each incarnation of that permutation (ie, only once for these
    % puppies)
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
    x = floor(log2(freq/200)/.1);
    y = floor(amp / 10);
    Z(x,y)= length(raster_arr);
end

figure('units','normalized','outerposition',[0 0 .25 .8]);

key = replace(key, '_', ' ');
title({key});%;['p2p: ' num2str(s.p2p) '            sym: ' num2str(s.sym)]})

surf(Z)
end
