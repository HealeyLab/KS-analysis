
function [fig, tosave,...
        BBavg, NNavg, BNavg, NBavg] = cross_corr_helper(obj,tet_cells,...
        BBavg, NNavg, BNavg, NBavg)
    db = obj.db;
    fig = figure('units','normalized','outerposition',[0 0 1 1]);
    hold on;

    tet_cells = obj.filter_tet_cells(tet_cells);
    tosave = 0;

    for i = 1:size(tet_cells,1)
        key_i = tet_cells(i,:);
        s_i = db(obj.keyhash(key_i{1}, key_i{2}, key_i{3}, key_i{4}));

        for j = i+1:size(tet_cells,1)

            key_j = tet_cells(j,:);
            s_j = db(obj.keyhash(key_j{1}, key_j{2}, key_j{3}, key_j{4}));

            [tsOffsets, ~, ~] = crosscorrelogram(...
                s_i.spike_timestamps / 1000, s_j.spike_timestamps / 1000, [-0.100 0.100]);

%                     subplot(length(tet_cells),length(tet_cells),...
%                         (i-1) * length(tet_cells) + j)   
%                     h = histogram(tsOffsets, 100);

            % set color                    
            BB = obj.BB; NN = obj.NN;
            BN = obj.BN; NB = obj.NB;
            set(gca, 'YTick', [], 'XTick', [])
            xlim([-.1, .1]);
            yl = ylim; line([0 0], [0 yl(2)], 'Color', 'k');

            tosave = isfield(s_i, 'p2p') && isfield(s_j, 'p2p');
            if tosave
                if s_i.p2p >= .43 && s_j.p2p >= .43
                    h.FaceColor = BB; h. EdgeColor = BB;
                elseif s_i.p2p < .43 && s_j.p2p < .43
                    h.FaceColor = NN; h. EdgeColor = NN;
                elseif s_i.p2p >= .43 && s_j.p2p < .43
                    h.FaceColor = BN; h. EdgeColor = BN;
                elseif s_i.p2p < .43 && s_j.p2p >= .43
                    h.FaceColor = NB; h.EdgeColor = NB;
                end

            %% quick analysis bit here at the end
            % will get devation from median. value of zero means
            % flat, value of negative meas inhibition, value of
            % positive means excitation

            % I only took from one side because though xcorr isn't
            % really commutative, you just don't see it playing out
            % that way in the data.
                if j > i
                    N = histcounts(tsOffsets, 100);
                    %TODO: MAKE THIS EVOKED ACTIVITY
                    % so this is median deviation. Wish me luck.
                    [lmi, l_min] = min(N(1:50));  [lma, l_max] = max(N(1:50)); 
                    [rmi, r_min] = min(N(51:100)); [rma, r_max] = max(N(51:100));
                    to_add = [l_min-50 l_max-50 r_min r_max]; % fix index

                    if s_i.p2p >= .43 && s_j.p2p >= .43
                        BBavg = [BBavg; to_add];
                    elseif s_i.p2p < .43 && s_j.p2p < .43
                        NNavg = [NNavg; to_add];
                    elseif s_i.p2p >= .43 && s_j.p2p < .43
%                                 if max_ind <= 2
                        BNavg = [BNavg; to_add];
%                                 else
%                                     NBavg = [NBavg; l_min l_max r_min r_max];
%                                 end
                    elseif s_i.p2p < .43 && s_j.p2p >= .43
%                                 if max_ind <= 2
                        NBavg = [NBavg; to_add];
%                                 else
%                                     BNavg = [BNavg; l_min l_max r_min r_max];
%                                 end
                    end
                end
            end
        end
    end      
hold off;
end
