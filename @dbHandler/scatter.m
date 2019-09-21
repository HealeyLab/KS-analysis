        function scatter(obj)
            vec=[];
            db = obj.db;
            for i = 1:length(obj.wf_keys)
                key = obj.wf_keys{i};
                s = db(key);
                disp(i)
                if strcmp(db(key).goodness, 'good')
                    if isfield(db(key), 'stim_timestamps')
                        [on, off, ~]=get_stim_on_off(obj, key);

                        total_spikes = 0;
                        total_time = 0;
                        for j = 2:length(db(key).stim_timestamps)-1 % ignore the first and last one because might be truncated.
                            num_sp = length(s.spike_timestamps(...
                                s.spike_timestamps>on(j) &...
                                s.spike_timestamps<off(j)));
                            total_spikes = total_spikes + num_sp;
                            total_time = total_time +...
                                1 / s.adc_sampling_rate * (off(j) - on(j));
                        end
                    else
                        spike_ts = db(key).spike_timestamps;
                        total_spikes = length(spike_ts);
                        total_time = (spike_ts(end) - spike_ts(1)) * 1 / s.amplifier_sampling_rate;
                    end
                    vec = [vec; s.p2p s.sym total_spikes/total_time s.depth];
                end
            end
            %scatter
            scatter(vec(:,1), vec(:,2), 'filled','k', 'jitter', 'on', 'jitterAmount', 0.01)
            xlabel('peak to peak width (ms)')
            ylabel('symmetry ratio')
            title('BS and NS neuron clustering')

            %scatter w depth
            figure;
            plotSpread(vec(:,4))
            scatter(vec(:,1), vec(:,4), 'filled','k', 'jitter', 'on', 'jitterAmount', 0.01)
            xlabel('peak to peak width (ms)')
            ylabel('depth (um)')
            title('width vs depth')
            
            %p2pvs fr
            figure;
            scatter(vec(:,1), vec(:,3), 'filled','k', 'jitter', 'on', 'jitterAmount', 0.01)
            xlabel('peak to peak width (ms)')
            ylabel('evoked firing rate (Hz)')
            title('Peak to peak vs firing rate')
            % k-means is garbage here
            idx = kmeans(vec(:,[1 3]), 3, 'Options', statset('Display','final'), 'Replicates', 5);
            scatter(vec(idx==1,1), vec(idx==1,3), 'MarkerEdgeColor', obj.BB, 'MarkerFaceColor',obj.BB, 'jitter', 'on', 'jitterAmount', 0.01)
            hold on;
            scatter(vec(idx==2,1), vec(idx==2,3), 'MarkerEdgeColor', obj.NN,'MarkerFaceColor', obj.NN, 'jitter', 'on', 'jitterAmount', 0.01)
%             scatter(vec(idx==3,1), vec(idx==3,3), 'MarkerEdgeColor', 'k','MarkerFaceColor', 'k')
%             
            %fr vs depth
            figure;
            scatter(vec(:,3), vec(:,4), 'filled','m', 'jitter', 'on')
            xlabel('Evoked firing rate (Hz)')
            ylabel('Depth (um)')
            title('fr vs depth')
            
            % show kmeans for p2p vs sym
            figure;
            idx = kmeans(vec(:,1:2), 2, 'MaxIter',100, 'Options', statset('Display','final'), 'Replicates', 5);
            scatter(vec(idx==1,1), vec(idx==1,2), 'MarkerEdgeColor', obj.BB, 'MarkerFaceColor',obj.BB)%, 'jitter', 'on')
            hold on;
            scatter(vec(idx==2,1), vec(idx==2,2), 'MarkerEdgeColor', obj.NN,'MarkerFaceColor', obj.NN)%, 'jitter', 'on')
            xlabel('peak to peak width (ms)')
            ylabel('symmetry ratio')
            line([0.43 0.43], [0 .7], 'Color', 'k');
            title('BS and NS neuron clustering')
            
            % p2p vs sym vs fr
            figure;
            scatter3(vec(:,1), vec(:,2), vec(:,3), 'jitter', 'on')
            xlabel('p2p')
            ylabel('symmetry')
            zlabel('fr')
            
        end
