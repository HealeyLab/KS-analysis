function [T, eva] = scatterTable(obj)

    T = table('Size', [1 7], 'VariableTypes',...
            {'string','double','double', 'double', 'double', 'double', 'double'},'VariableNames',...
            {'key'   ,'p2p',   'sym',    'FR',     'depth',  'TOD',    'id'});
%               1      2         3          4         5         6        7 

    db = obj.db;
    keys = obj.get_keys('&', 'trial'); % using & as a pattern excludes the auxiliary entries
    for key_ind = 1:length(keys)
        key = keys{key_ind};
        s = db(key);

        % filtering for ones that are good and for wf analysis
        if isfield(s, 'wf_analysis_include') && s.wf_analysis_include
            % Get evoked firing rate, if applicable
            total_spikes = 0;
            total_time = 0;    
            if isfield(db(key), 'stim_timestamps')
                [on, off, ~]=get_stim_on_off(obj, key);

                for j = 2:length(db(key).stim_timestamps)-1 % ignore the first and last one because might be truncated.
                    num_sp = length(s.spike_timestamps(...
                        s.spike_timestamps>on(j) &...
                        s.spike_timestamps<off(j)));
                    total_spikes = total_spikes + num_sp;
                    total_time = total_time +...
                        1 / s.adc_sampling_rate * (off(j) - on(j));
                end
            else
                total_spikes = NaN; total_time = NaN;
            end
            [~, id_num] = obj.get_subject_id(key);
            
            
            T = [T; {key obj.get_p2p(s) obj.get_sym(s) total_spikes/total_time s.depth   obj.get_time_of_day(key) id_num}];

            disp([num2str(height(T)) ' ' key])    
        end
    end
    T(1,:)=[]; % removes first crap entry
    % Finally, cluster them and add that to the tbale
    eva = evalclusters([T.p2p T.sym],'kmeans', 'Gap', 'KList',[1:5]);
    
    Z = linkage([rescale(T.p2p) rescale(T.sym)], 'Ward');
    hclust = cluster(Z, 'maxclust' ,eva.OptimalK);
    km = kmeans([rescale(T.p2p) rescale(T.sym)], eva.OptimalK, 'MaxIter',100, 'Options', statset('Display','final'), 'Replicates', 5);
    t = table(hclust, km);
    T = [T t ];
    

end