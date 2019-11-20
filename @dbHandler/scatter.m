function vec = scatter(obj)
    % This function deals a lot with whether you are using a particular
    % unit for waveform analysis. Let's get into what that means.
    % for_wf_analysis: is for this kind of analysis, w/o narrowing down
    % included: final pass. yes or no, this is going to be used in
    % the analysis.
    % goodness: 'good ', 'mua  ', or 'noise'. Note the spaces. How it was
    % sorted in Kilosort.
    if ~exist('vec', 'var')
        vec=[];
        db = obj.db;
        keys = obj.get_keys('&', 'trial'); % using & as a pattern excludes the auxiliary entries
        for key_ind = 1:length(keys)
            key = keys{key_ind};
            s = db(key);

            % filtering for ones that are good and for wf analysis
            if strcmp(strrep(db(key).goodness,' ', ''), 'good') && db(key).for_wf_analysis

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
                
                %                   1              2              3                  4                     5           6
                vec = [vec; obj.get_p2p(s) obj.get_sym(s) total_spikes/total_time s.depth   obj.get_time_of_day(key) id_num];
                
                disp([num2str(length(vec)) ' ' key])    
            end
            

        end
    end
    %% show kmeans for p2p vs sym
    subplot(3,2,1)
    idx = kmeans(vec(:,1:2), 2, 'MaxIter',100, 'Options', statset('Display','final'), 'Replicates', 5);
    scatter(vec(idx==1,1), vec(idx==1,2), 'MarkerEdgeColor', obj.BB, 'MarkerFaceColor',obj.BB)%, 'jitter', 'on')
    hold on;
    scatter(vec(idx==2,1), vec(idx==2,2), 'MarkerEdgeColor', obj.NN,'MarkerFaceColor', obj.NN)%, 'jitter', 'on')
    xlabel('peak to peak width (ms)')
    ylabel('symmetry ratio')
    line([0.43 0.43], [0 .7], 'Color', 'k');
    title('BS and NS neuron clustering')
    
    %% scatter w depth
    subplot(3,2,2)
%     plotSpread(vec(:,4))
    scatter(vec(:,1), vec(:,4), 'filled','k', 'jitter', 'on', 'jitterAmount', 0.01)
    xlabel('peak to peak width (ms)')
    ylabel('depth (um)')
    title('width vs depth')

    %% p2p vs fr
    subplot(3,2,3)
    scatter(vec(:,1), vec(:,3), 'filled','k', 'jitter', 'on', 'jitterAmount', 0.01)
    xlabel('peak to peak width (ms)')
    ylabel('evoked firing rate (Hz)')
    title('Peak to peak vs evoked firing rate')
    %% depth vs identity
    subplot(3,2,4)
    scatcell = cell(1,1);
    for i = 1:8
        scatcell{i} = vec(vec(:,6)==i,4); % depth is 4, id is 6
    end
    catLabels = {'mda','mdb','mdc','mdd','mde','mdx','mdy','mdz'};
    plotSpread(scatcell, 'CategoryLabels', catLabels)
    %     scatter(vec(:,3), vec(:,4), 'filled','m', 'jitter', 'on')
    xlabel('Evoked firing rate (Hz)')
    ylabel('Depth (um)')
    title('Depth vs identity')

    %% fr vs time of day
    subplot(3,2,5)
    height = max(vec(:,3)+10);
    xlim([0 24])
    ylim([0 height])
    hold on
    area([0:6],ones(7,1) * height, 'FaceColor', [0.8 0.8 0.8], 'LineStyle','None');
    area([18:24],ones(7,1) * height, 'FaceColor', [0.8 0.8 0.8], 'LineStyle','None');
    scatter(vec(:,5), vec(:,3), 'filled', 'k', 'jitter', 'on');
    title('FR vs time of day');
    xlabel('time (hr)')
    ylabel('Evoked Fr (Hz)')    
    %% stuff vs hemisphere
    
    %     %% p2p vs sym vs fr
%     subplot(3,2,6)
%     scatter3(vec(:,1), vec(:,2), vec(:,3), 'jitter', 'on')
%     title('p2p vs symmetry vs evoked fr')
%     xlabel('p2p')
%     ylabel('symmetry')
%     zlabel('evoked fr (Hz)')
end
