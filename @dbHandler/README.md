#### dbHandler.add
Adds all the units in a recording to the database. Overwrites existing key entries FOR ALL UNITS IN A RECORDING if key already exists, i.e., if you've added that recording to the database before and you wanted to change its entry.

#### dbHandler.extract_stim_timestamps
Take precise timestamps for stimulus onsets to be aligned with stimulus waveforms later. Uses cross-correlation to find the exact onset of a .wav file in the recording.

#### dbHandler.clean_peaks
Internal function to aid `dbh.extract_stim_timestamps`. Prevents redundancies in thresholding onsets of stimuli.

#### dbHandler.filter_raw_data
The current version of the KiloSort spike-sorting script filters and saves your data for you, so this should be unnecessary for most applications. NOTE: Make sure the current directory contains the raw .rhd data before running this function.

#### dbHandler.getWaveFormsDriver
Uses a the cortex-lab's Spikes repository (which should be installed) to grab waveforms from your spike data. Make sure you are in the Kilosort folder when you run this.

#### dbHandler.get_keys(obj, pattern, exclude)
Returns the keys that match the pattern (all if no parameters given)

#### dbHandler.get_key_family(obj, key)
Keys are in the form: `'mde 10 13 19 fem return_191013_172903_Kilosort&96&13&good'`. This function strips the identifying characteristics of the specific unit, yielding `'mde 10 13 19 fem return_191013_172903_Kilosort'`. It then searches the other keys in the database for matches, and returns them as a cell. Note: recent analyses have created necessity for non-typical keys. 

#### dbHandler.gen_playback_syllable_PSTHs(obj, key)
Uses syllable onsets and offsets you obtained using `dbHandler.get_playback_syllables(obj)` to generate a PSTH for one cell for each syllable in each song for each presentation of that song. It also returns the timestamps plotted in the PSTH.

#### dbHandler.get_playback_syllables(obj, key)
So you've done your recording. You've sorted your data. Now you want to know how your cells respond to syllables in your playback files. Manually hardcode the `subject` and `path` variables in the function to point to your zf son folder, and then use the gui to pick out each syllable. For BOS and BOS-REV, you may want to take care to capture exactly the same syllables.

#### dbHandler.get_song_syllable_activity(obj, key, only)
Generates syllable PSTHs based on user input to a gui. The `only` parameter is for if you want to generate PSTHs for all cells in the key's "family" (see dbHandler.get_key_family)

#### dbHandler.keyhash(~, workingDirectory, unit, channel, goodness)
Produces a key for the containers.Map object to use for indexing a particular unit. The key includes the name of the recording, the time of the recording, the channel, and the cluster. Feel free to add more information to the key if your procedure demands it.

#### dbHandler.plot_waveform(obj,key, varargin)
Plots a waveform on arbitrary given axis (varargin argument) as broad or narrow (based on our arbitrary threshold), including grey background of 100 individual waveforms, the average, plotted in the color of broad of narrow cells, and plus or minus stdev.

#### dbHandler.scatter.m
Gathers data for and generates summary figures

#### dbhHandler.waveform_analysis(obj)
Generates a gui for making the last pass for waveform analysis. Users label each waveforms as "include" or "exclude"

#### dbHandler.waveform_connector(obj,habit_key, song_key)
For the "key families" of both parameters, generates a subplot graph of each cell. The user must visually inspect each one and determine which waveforms correspond, taking into account whether the waveforms look the same and are on the same tetrode. They then input the indicies (indicated in each subplot) of each pair of waveforms they believe to correspond, one at a time. When all waveforms have been matched, enter "done". A gui will pop up, prompting you to annotate a vocalization. Once that is done, the function will save the pairings and timestamps in a file named latency_db.mat. To access the database, I recommend the pattern:

```
latency_db = load('PATH\TO\LATENCY\DB.MAT');
latency_db = latency_db.latency_db
```

This is a heirarchical database implemented naively with structs. No fancy low-memory use table implementations here. The structure is:
latency_db(cell1(song1(sylA({timestamps}) sylB({timestamps}) ...)...) cell2(...)...)
Code for summarizing this database's contents is forthcoming.

