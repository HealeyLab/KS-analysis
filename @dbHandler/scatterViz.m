function scatterViz(obj, T,eva)
%% SCATTER 
%% show kmeans for p2p vs sym
colors = {obj.BB obj.BN obj.NB obj.NN};
clust = {'hclust','km'};
for i = 1:2
figure('Name',clust{i})
    for ii = 1:length(unique(T.(clust{i}))) % numof clusts
        color = colors{ii};
        
        scatter(T.p2p(T.(clust{i})==ii), T.sym(T.(clust{i})==ii), 'MarkerEdgeColor', color,...
            'MarkerFaceColor',color) %, 'jitter', 'on','jitterAmount', 0.01)
        hold on;
    end
end
xlabel('peak to peak width (ms)')
ylabel('symmetry ratio')
line([obj.p2p_thres obj.p2p_thres], [0 .7], 'Color', 'k');
title('BS and NS neuron clustering')

subplot(3,2,1)
gscatter(T.p2p, T.sym, eva.OptimalY, 'rbgkc','xod^*')

%% scatter w depth
subplot(3,2,2)
scatter(T.p2p, T.depth, 'filled','k')%, 'jitter', 'on', 'jitterAmount', 0.01)
xlabel('peak to peak width (ms)')
ylabel('depth (um)')
title('width vs depth')

%% p2p vs fr
subplot(3,2,3)
scatter(T.p2p, T.FR, 'filled','k', 'jitter', 'on', 'jitterAmount', 0.01)
xlabel('peak to peak width (ms)')
ylabel('evoked firing rate (Hz)')
title('Peak to peak vs evoked firing rate')
%% depth vs identity
subplot(3,2,4)
scatcell = cell(1,1);
for i = 1:8
    scatcell{i} = T.depth(T.id==i);
end
catLabels = {'mda','mdb','mdc','mdd','mde','mdx','mdy','mdz'};
plotSpread(scatcell, 'CategoryLabels', catLabels)
xlabel('Subject #')
ylabel('Depth (um)')
title('Depth vs identity')

%% fr vs time of day
subplot(3,2,5)
height = max(T.FR)+10;
xlim([0 24])
ylim([0 height])
hold on
area([0:6],ones(7,1) * height, 'FaceColor', [0.8 0.8 0.8], 'LineStyle','None');
area([18:24],ones(7,1) * height, 'FaceColor', [0.8 0.8 0.8], 'LineStyle','None');
scatter(T.TOD, T.FR, 'filled', 'k', 'jitter', 'on');
title('FR vs time of day');
xlabel('time (hr)')
ylabel('Evoked Fr (Hz)')    
%% stuff vs hemisphere

%% p2p vs sym vs fr
subplot(3,2,6)
scatter3(T.p2p,T.sym, T.FR, 'jitter', 'on')
title('p2p vs symmetry vs evoked fr')
xlabel('p2p')
ylabel('symmetry')
zlabel('evoked fr (Hz)')
end
