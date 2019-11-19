function [outputArg1,outputArg2] = plot_waveform(obj,key, varargin)
%PLOT_WAVEFORM I just modularized the plotting of waveforms
%   Now, this can be called on an arbitrary axis
if ~isempty(varargin)
    ax = varargin{1};
    set(gcf,'CurrentAxes',ax); 
end

s = obj.db(key);
wf = s.spike_waveforms;
wf_mean = nanmean(wf,2);
hold on

s.p2p = obj.get_p2p(s);
s.sym = obj.get_sym(s);

plot(wf(:,1:100), 'y', 'Color', [0,0,0]+.75)
wf_std = nanstd(wf')';
plot(wf_mean, 'LineWidth', 1.25, 'Color', obj.get_color(s));
plot(wf_mean+wf_std,'LineWidth', 1.25, 'Color', 'k');
plot(wf_mean-wf_std,'LineWidth', 1.25, 'Color', 'k');  
set(gca, 'XTick', [])
ylim([(min(wf_mean-wf_std*5)) (max(wf_mean+wf_std)*5)]);

end

