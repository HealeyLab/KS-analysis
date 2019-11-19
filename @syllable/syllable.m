classdef syllable
    %SYLLABLE Holds data about waveform and neural activity of a syllable
    %   Stores two pieces of information: 
    %   window_s    : timestamps at beginning and end of sound
    %   cells       : a containers.Map that contains a struct with the following:
    %                    spiketrain     : the spikes during the sound
    %                    key            : keyhash for that cell
    
    properties (SetAccess = immutable)
        window_s;
        cells     = {};
        cell_keys = {}; % for cross referencing the order
        id = ''
        sonogram;
        adc_sr;
        amplifier_sr;
    end
    
    methods
        function obj = syllable(dbh, window_s, cell_keys, id, varargin)
            %SYLLABLE Construct an instance of this class
            % keys is a cell array of keyhashes    
            obj.window_s = window_s; % x(1) and x(2)
            obj.id = id;
            obj.cell_keys = cell_keys;
            % uses the keyhash to grab the spiketrain, and narrows down
            % only the spikes that happen during the syllable, using the
            % window_s's timestamps.
            
            %% 
            % first, convert window_s's timestamps to amplifier timestamps
            for i=1:length(cell_keys)
                % get entry and find timestamps
                entry = dbh.db(cell_keys{i});
                
                if i == 1
                    % sono is the syllable.window_s field, a [double double]
                    mic = dbh.get_microphone(cell_keys{1});
                    obj.adc_sr = entry.adc_sampling_rate;
                    obj.amplifier_sr = entry.amplifier_sampling_rate;
                    win = uint64(obj.window_s * obj.adc_sr * 60); % convert back to samples
                    sonogram = mic(win(1):win(2));
                    obj.sonogram = sonogram - mean(sonogram);
                end
                
                amp_x = obj.window_s * 60 * entry.amplifier_sampling_rate;
                sTs = entry.spike_timestamps; % query value
                
                % only get within this window and normalize
                obj.cells{i} = ...
                    sTs(sTs > amp_x(1) & sTs < amp_x(2)) - amp_x(1); 
                
                if exist(varargin, 'var')
                    obj.sonogram = varargin{1};
                end
            end            
        end
    end
end

