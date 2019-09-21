function waveform_analysis(obj)
    fig = figure('Visible','on','Position',[360,200,550,300]);

    values = uicontrol('Style', 'text', 'String', '',...
                'Tag', 'p2psym',...
                'Position', [130, 10, 140, 25]);

    type = uicontrol('Style', 'text', 'String', '',...
                'Tag', 'type',...
                'Position', [315, 240, 70, 25]);

    status = uicontrol('Style', 'text', 'String', '',...
                'Tag', 'status',...
                'Position', [315, 200, 70, 25]);

    progress = uicontrol('Style', 'text', 'String', '',...
                'Tag', 'progress',...
                'Position', [315, 170, 70, 25]);

    good    = uicontrol('Style','pushbutton',...
                 'String','Good',...
                 'Position',[315,140,70,25],...
                 'Callback',@good_Callback);

    bad    = uicontrol('Style','pushbutton',...
                 'String','Bad','Position',[315,100,70,25],...
                 'Callback',@bad_Callback);

    next = uicontrol('Style','pushbutton',...
                 'String','next','Position',[315,60,70,25],...
                 'Callback', @next_Callback);

    back = uicontrol('Style','pushbutton',...
                 'String','back','Position',[315,20,70,25],...
                 'Callback', @back_Callback);

    axes('Units','Pixels','Position',[50,60,200,185]); 

    psth = uicontrol('Style','pushbutton',...
                    'String', 'PSTH',...
                    'Position', [400, 60, 70, 25],...
                    'Callback', @psth_Callback);

    xcorr = uicontrol('Style', 'pushbutton',...
        'String', 'xcorr',...
        'Position', [400, 100, 70, 25],...
        'Callback', @xcorr_Callback);

    function good_Callback(source,~) 
    % for if is good
        [s,key] = get_s;
        s.wf_analysis_goodness = 'Good';
        obj.db(key) = s;
        handle = ancestor(source, 'figure');
        status_handle = findobj(handle, 'Tag', 'status');
        set(status_handle, 'String', 'Included: good');                
    end

    function bad_Callback(source,~) 
    % for if is bad
        [s, key] = get_s;
        s.wf_analysis_goodness = 'Bad';
        obj.db(key) = s;
        handle = ancestor(source, 'figure');
        status_handle = findobj(handle, 'Tag', 'status');
        set(status_handle, 'String', 'Included: bad');
    end

    function back_Callback(source,~) 
    % Display contour plot of the currently selected data.
        handle = ancestor(source, 'figure');
        s = increment(-1, handle);

        status_handle = findobj(handle, 'Tag', 'status');
        if isfield(s, 'wf_analysis_goodness')
            set(status_handle, 'String', s.wf_analysis_goodness);
        else
            set(status_handle, 'String', 'Included: undecided');
        end 
    end

    function next_Callback(source,~) 
        % for getting next
        handle = ancestor(source, 'figure');

        s = increment(1, handle);

        status_handle = findobj(handle, 'Tag', 'status');
        if isfield(s, 'wf_analysis_goodness')
            set(status_handle, 'String', s.wf_analysis_goodness);
        else
            set(status_handle, 'String', 'Included: undecided');
        end        
    end        

    function s = increment(amt, handle)
        obj.count = obj.count+amt; 

        [s, key] = get_s;

        wf = s.spike_waveforms;
        wf_mean = nanmean(wf,2);
        plot(wf_mean, 'LineWidth', 2); hold on;
        plot(wf_mean+nanstd(wf')', 'b');
        plot(wf_mean-nanstd(wf')', 'b');  

        type_handle = findobj(handle, 'Tag', 'type');
        set(type_handle, 'String', ['Sorted as:' s.goodness]);

        type_handle = findobj(handle, 'Tag', 'progress');
        set(type_handle, 'String', [num2str(obj.count) '/' num2str(length(obj.wf_keys))]);

        [Vmax,Imax] = max(wf_mean);
        [Vmin,Imin] = min(wf_mean);
        s.p2p = abs(Imax - Imin) / s.amplifier_sampling_rate * 1000 ;
        s.sym = abs(Vmax/Vmin);
        plot(Imax, Vmax, 'm*');
        plot(Imin, Vmin, 'm*');
        hold off;
        obj.db(key) = s;

        values_handle = findobj(handle, 'Tag', 'p2psym');
        set(values_handle, 'String', ['p2p: ' num2str(s.p2p) '  sym: ' num2str(s.sym)]);
    end

    function [s,key] = get_s
        keys = obj.wf_keys;
        key = keys{obj.count};
        s = obj.db(key);
    end


    function psth_Callback(~,~)
        % will eventually be the callback for the save function
        db = obj.db; keys = obj.wf_keys;
        key = keys{obj.count};
        s = db(key);
        [on, off, wav_files] = obj.get_stim_on_off(key);

        % for each class of stim:
        si = s.stim_identities{1};
        usi = unique(si);
        for i = 1:length(usi)
            % for each stim itself
            stim_inds = find(strcmp(si,usi{i}));
            raster_arr = cell(length(stim_inds),1);
            for j = 1:length(stim_inds)
                sp_ts = s.spike_timestamps;
                ind = stim_inds(j);
                start = on(ind) - 2 * s.amplifier_sampling_rate;
                stop = off(ind) + 2 * s.amplifier_sampling_rate;

                raster_arr{j} = ((intersect(...
                    sp_ts(sp_ts > start),...
                    sp_ts(sp_ts < stop))...
                    - start) / s.amplifier_sampling_rate)';                        
            end

            figure;
            subplot(3,1,1);
            for k = 1:length(wav_files)
                if contains(usi{i}, wav_files(k).name)
                    curr_wav = wav_files(k);
                end
            end

            spectrogram(... 
                [nan(2*s.adc_sampling_rate, 1); curr_wav.data; nan(2*s.adc_sampling_rate, 1)],...
                256, [],[],s.adc_sampling_rate, 'yaxis')
            colorbar('delete');

            title([curr_wav.name '    ' s.context '  why on earth does this histogram work']);

            subplot(3,1,2);
            [xpoints, ~] = plotSpikeRaster(raster_arr,...
                    'PlotTYpe','vertline', 'XLimForCell', [0 (stop-start)/s.amplifier_sampling_rate]);

            histo_axes = subplot(3,1,3);    
            bin = 0.010; % bin size in s
            histogram(histo_axes, xpoints, (0:bin:(stop-start)/s.amplifier_sampling_rate)); % convert ms to s
            xlim([0, (stop-start)/s.amplifier_sampling_rate])
            % construct psth for that stim
        end
    end

    function xcorr_Callback(~,~)
        [~,key] = get_s;
        nothin = obj.cross_correlogram_for_gui(key);
    end
end
