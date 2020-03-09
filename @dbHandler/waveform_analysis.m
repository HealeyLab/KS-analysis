function waveform_analysis(obj)
%%WAVEFORM_ANALYSIS This is for going through your entire dataset and 
%%marking each waveform as good or bad. It is your final pass.
% This opens a GUI to let you increment through each waveform. It displays
% some statistics about the waveform, including p2p length and symmtry.
    keys = obj.db.keys;
    wf_keys = {};
    for i = 1:length(keys)
        s = obj.db(keys{i});
        if isfield(s, 'for_wf_analysis') && s.for_wf_analysis
            s.wf_analysis_include = 0;
            wf_keys{length(wf_keys) + 1} = keys{i};
        end
    end
    count = 0;
    
    fig = figure('Visible','on','Position',[360,200,550,300]);

    values = uicontrol('Style', 'text', 'String', '',...
                'Tag', 'p2psym','Position', [130, 10, 140, 25]);

    type = uicontrol('Style', 'text', 'String', '',...
                'Tag', 'type','Position', [315, 240, 70, 25]);

    status = uicontrol('Style', 'text', 'String', '',...
                'Tag', 'status',...
                'Position', [315, 200, 70, 25]);

    progress = uicontrol('Style', 'text', 'String', '',...
                'Tag', 'progress','Position', [315, 170, 70, 25]);

    good    = uicontrol('Style','pushbutton', 'String','Include',...
                 'Position',[315,140,70,25], 'Callback',@include_Callback);

    bad    = uicontrol('Style','pushbutton','String','Exclude',...
                 'Position',[315,100,70,25], 'Callback',@exclude_Callback);

    next = uicontrol('Style','pushbutton','String','next',...
                 'Position',[315,60,70,25],'Callback', @next_Callback);

    back = uicontrol('Style','pushbutton','String','back',...
                 'Position',[315,20,70,25], 'Callback', @back_Callback);

    axes('Units','Pixels','Position',[50,60,200,185]); 

    psth = uicontrol('Style','pushbutton', 'String', 'PSTH',...
                    'Position', [400, 60, 70, 25],'Callback', @psth_Callback);

    function include_Callback(source,~) 
    % for if is good
        [s,key] = get_s();
        s.wf_analysis_include = 1;
        obj.db(key) = s;
        handle = ancestor(source, 'figure');
        status_handle = findobj(handle, 'Tag', 'status');
        set(status_handle, 'String', 'INCLUDED');                
    end

    function exclude_Callback(source,~) 
    % for if is bad
        [s, key] = get_s();
        s.wf_analysis_include = 0;
        obj.db(key) = s;
        handle = ancestor(source, 'figure');
        status_handle = findobj(handle, 'Tag', 'status');
        set(status_handle, 'String', 'UNINCLUDED');
    end

    function show_wf_analysis_include(s, status_handle)
        if isfield(s, 'wf_analysis_include') && s.wf_analysis_include
            include = 'INCLUDED';
        else
            include = 'NOT INCLUDED';
        end
        
        if isfield(s, 'wf_analysis_include')
            set(status_handle, 'String', include);
        else
            set(status_handle, 'String', 'undecided');
        end 
    end

    function back_Callback(source,~) 
    % Display contour plot of the currently selected data.
        handle = ancestor(source, 'figure');
        s = increment(-1, handle);
        status_handle = findobj(handle, 'Tag', 'status');
        show_wf_analysis_include(s, status_handle);
    end

    function next_Callback(source,~) 
        % for getting next
        handle = ancestor(source, 'figure');
        s = increment(1, handle);
        status_handle = findobj(handle, 'Tag', 'status');
        show_wf_analysis_include(s, status_handle);
    end        

    function s = increment(amt, handle)
        count = count+amt; 

        [s, key] = get_s();

        wf = s.spike_waveforms;
        wf_mean = nanmean(wf,2);
        % plot 20 in grey
        plot(wf(:, randperm(2000,min(50,size(wf,2)))), 'Color', [.7 .7 .7]); hold on;
        % plot the mean
        plot(wf_mean, 'LineWidth', 2);
        % plot the std plus or minus
        plot(wf_mean+nanstd(wf')', 'b');
        plot(wf_mean-nanstd(wf')', 'b');  

        type_handle = findobj(handle, 'Tag', 'type');
        set(type_handle, 'String', ['Sorted as:' s.goodness]);

        type_handle = findobj(handle, 'Tag', 'progress');
        set(type_handle, 'String', [num2str(count) '/' num2str(length(wf_keys))]);

        
        [Vmin,Imin] = min(wf_mean);
        [Vmax,Imax] = max(wf_mean(Imin:end));
        s.p2p = obj.get_p2p(s);
        s.sym = obj.get_sym(s); % abs(Vmax/Vmin);
        plot(Imax+Imin, Vmax, 'm*');
        plot(Imin, Vmin, 'm*');
        ylim([Vmin*3 Vmax*7])
        hold off;
        obj.db(key) = s;

        values_handle = findobj(handle, 'Tag', 'p2psym');
        set(values_handle, 'String', ['p2p: ' num2str(s.p2p) '  sym: ' num2str(s.sym)]);
    end

    function [s,key] = get_s()
        key = wf_keys{count};
        s = obj.db(key);
    end

    function psth_Callback(~,~)
        % will eventually be the callback for the save function
        key = wf_keys{count};
        obj.gen_psth(key)
    end
end
