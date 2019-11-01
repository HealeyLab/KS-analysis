classdef syllable
    %SYLLABLE Holds data about waveform and neural activity of a syllable
    %   Stores two pieces of information: 
    %   sonogram    : timestamps at beginning and end of sound
    %   cells       : a containers.Map that contains a struct with the following:
    %                    spiketrain     : the spikes during the sound
    %                    key            : keyhash for that cell
    
    properties
        dbh;
        sonogram;
        cells = containers.Map('KeyType','char','ValueType', 'any');
        id = ''
    end
    
    methods
        function obj = syllable(dbh, sonogram, keys, id)
            %SYLLABLE Construct an instance of this class
            % keys is a cell array of keyhashes    
            obj.dbh = dbh;
            obj.sonogram = sonogram; % x(1) and x(2)
            obj.id = id;
            % uses the keyhash to grab the spiketrain, and narrows down
            % only the spikes that happen during the syllable, using the
            % sonogram's timestamps.
            
            % first, convert sonogram's timestamps to amplifier timestamps
            for i=1:length(keys)
                entry = dbh.db(keys{i});
                amp_x = obj.sonogram * 60 * entry.amplifier_sampling_rate;
                sTs = entry.spike_timestamps; % query value
                obj.cells(keys{i}) = sTs(sTs > amp_x(1) & sTs < amp_x(2)); % only get within this window
            end
        end
    end
end

