function [fig] = cross_corr_helper(obj,tet_cells)
    db = obj.db;
    fig = figure('units','normalized','outerposition',[0 0 1 1]);
    hold on;
    tet_cells = obj.filter_tet_cells(tet_cells);
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

            p2p_i = obj.get_p2p(s_i);
            p2p_j = obj.get_p2p(s_j);

            if     p2p_i >= .43&& p2p_j >= .43
                h.FaceColor = BB; h. EdgeColor = BB;
            elseif p2p_i < .43  && p2p_j < .43
                h.FaceColor = NN; h. EdgeColor = NN;
            elseif p2p_i >= .43 && p2p_j < .43
                h.FaceColor = BN; h. EdgeColor = BN;
            elseif p2p_i < .43  && p2p_j >= .43
                h.FaceColor = NB; h.EdgeColor = NB;
            end

        end
    end      
hold off;
end
