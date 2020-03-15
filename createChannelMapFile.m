%  create a channel map file

Nchannels = 32;
connected = true(Nchannels, 1);
chanMap   = 1:Nchannels;

chanMap0ind = chanMap - 1;
xcoords   = ones(Nchannels,1);
ycoords   = [1:Nchannels]';
kcoords   = ones(Nchannels,1); % grouping of channels (i.e. tetrode groups)
fs = 25000; % sampling frequency

% the first thing Kilosort does is reorder the data with data = data(chanMap, :).
% Now we declare which channels are "connected" in this normal ordering, 
% meaning not dead or used for non-ephys data

connected = logical([zeros(8,1); ones(16,1); zeros(8,1)]);
% now we define the horizontal (x) and vertical (y) coordinates of these
% 34 channels. For dead or nonephys channels the values won't matter. Again
% I will take this information from the specifications of the probe. These
% are in um here, but the absolute scaling doesn't really matter in the
% algorithm. 

xcoords = repmat([1 2 3 4]', 1, Nchannels/4);
xcoords = xcoords(:);

% Often, multi-shank probes or tetrodes will be organized into groups of
% channels that cannot possibly share spikes with the rest of the probe. This helps
% the algorithm discard noisy templates shared across groups. In
% this case, we set kcoords to indicate which group the channel belongs to.
% In our case all channels are on the same shank in a single group so we
% assign them all to group 1. 

number_of_tetrodes = 4;
tetrode1 = [8:11]; tetrode2 = [12:15];
tetrode3 = [16:19]; tetrode4 = [20:23];

all_tetrodes = [tetrode1;tetrode2;tetrode3;tetrode4] + 1;

kcoords =  ones(Nchannels, 1);
for x=1:number_of_tetrodes
    kcoords(all_tetrodes(x,:)) = x;
    xcoords(all_tetrodes(x,:)) = [1:4];
end
xcoords = xcoords(:);
ycoords = kcoords;

save('C:\DATA\Spikes\20150601_chan32_4_900s\chanMap.mat', ...
    'chanMap','connected', 'xcoords', 'ycoords', 'kcoords', 'chanMap0ind', 'fs')
%%

% kcoords is used to forcefully restrict templates to channels in the same
% channel group. An option can be set in the master_file to allow a fraction 
% of all templates to span more channel groups, so that they can capture shared 
% noise across all channels. This option is

% ops.criterionNoiseChannels = 0.2; 

% if this number is less than 1, it will be treated as a fraction of the total number of clusters

% if this number is larger than 1, it will be treated as the "effective
% number" of channel groups at which to set the threshold. So if a template
% occupies more than this many channel groups, it will not be restricted to
% a single channel group. 