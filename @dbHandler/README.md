### dbh.add
Adds all the units in a recording to the database. Overwrites existing key entries FOR ALL UNITS IN A RECORDING if key already exists, i.e., if you've added that recording to the database before and you wanted to change its entry.

### dbh.extract_stim_timestamps
Like it says in the name, it gets precise timestamps for stimulus onsets. Can be aligned with stimulus waveforms later.

### dbh.clean_peaks
Internal function to aid `dbh.extract_stim_timestamps`. Prevents redundancies in thresholding.

### dbh.filter_raw_data
I made this because I had to filter a lot of my data post-hoc. The current version of the script filters and saves your data for you, so this should be unnecessary for most applications. YOu should be in the folder containing the raw .rhd data when you run this.

### dbh.getWaveFormsDriver
Uses a the cortex-lab's spikes repository, which you should have installed, to grab waveforms from your spike data. Make sure you are in the Kilosort folder when you run this.

### dbh.keyhash
Produces a key for the containers.Map object to use for indexing a particular unit. The key includes the name of the recording, the time of the recording, the channel, and the cluster. Feel free to add more information to the key if your procedure demands it.
