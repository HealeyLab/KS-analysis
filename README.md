# KS-analysis
A MATLAB and Python framework for unpacking KiloSort/phy data, putting it into a database, and visualizing it.

This small repo is intended to be used in conjunction with zeebie15's [ephysSuite](https://github.com/zeebie15/ephysSuite) repo, since it makes some assumptions about how your data is structured. See that repository for instructions for how to set up your rig to output the right data format.

## Installation
This comes with a few dependencies. Assuming you've successfully installed KiloSort and Phy, you'll also need to install the MATLAB engine for python. For that, you'll need another Python installation alongside the miniconda/Anaconda one you've got from Phy/KiloSort: 

https://www.python.org/downloads/release/python-369/

Follow the default install configuration. Then, you can invoke python 3.6 (the engine only supports 3.6, 3.5, and 2.7) while running your template-gui from 3.7+! To do so, open Anaconda prompt with Administrator priviledges and run `py -3.6 setup.py install`

https://www.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html

By the way, to install other dependencies in this console, just type 

``py -3.6 -m pip install <package>``

## Running Kilosort
Kilosort seems to demand a lot of manual copy/pasting of entire files, and editing of parameters and such by hand. I tried to automate most of this in one script, which should take minimal editing to get running on your machine. Yes, this will involve hard-coding. It might occur to you that you could just hardcode the files directly instead of through this script. I prefer doing things through the script because it allows me to leverage the (more or less) original templates kilosort provides from one place, without changing the default shape of the template files.

First, you'll have to edit the line reading
```Matlab
working_dir = 'C:\Users\danpo\Documents\MATLAB\DJP_KiloSort';
```
to 

```Matlab
working_dir = 'path\to\this\repos\directory';
```

Next, you'll have to add two lines to the part of the code that edits the master_file_example_MOVEME that point to your installation of Kilosort and to your installation of npy-matlab.

```Matlab
C{6} = sprintf('pathToYourConfigFile = ''%s''; ', dataPath);
% add the follwing lines below, with the indicated paths hardcoded:
C{3} = addpath(genpath('PATH\TO\YOUR\KILOSORT\INSTALL')) % path to kilosort folder
C{4} = addpath(genpath('PATH\TO\YOUR\NPY-MATLAB\INSTALL')) % path to npy-matlab scripts
```

The defaults I have set in `params.py`, `master_file_example_MOVEME.m`, `StandardConfig_MOVEME.m`, and `createChannelMapFile.m` are for a linear array of four tetrodes. On the Intan RHD2000 system, which records 32 channels, I only use channels 9 through 24 (16 total), and these files specify this configuration. I encourage you to look at them, try to get an intuition for how they work, and adjust them to your electrode configuration.

To run the script, make sure you have Kilosort and phy installed and are in the same directory as your .rhd files and run `DJP_KiloSort_Driver.m`.

Once you've run Kilosort, make sure to add the right \*markers.txt and \*stimtimes.txt files to the Kilosort folder. To check for the right pair of folders, make sure the first and last timestamp in the \*stimtimes.txt file corresponds to the beginning and end of the recording.

## Running analysis
To analyze your data, I designed a class named dbHandler for manipulating a MATLAB container.Map object that serves as your databse. Before you can use this class, you have to save a .mat file containing an empty container.Map object by running the following code in the directory you want to save your database:

```Matlab
db = containers.Map('KeyType','char','ValueType', 'any');
save('db.mat', 'db','-v7.3')
```

Make sure you have the [spikes repository](https://github.com/cortex-lab/spikes/) installed

I decided to hardcode in several useful paths into my class to make it easier to initialize the handler object across MATLAB sessions. So, change the properties section of the class to fit your needs:

```Matlab
properties
    dbPath = 'C:\Users\danpo\Documents\MATLAB\db.mat'; % --> C:\path\to\your\db.mat
    
    % Two different audio paths for two different experimental subjects, mdx and mdy
    audioPathMDX = 'C:\Users\danpo\Documents\MATLAB\ephysSuite\zf son mdx';
    audioPathMDY = 'C:\Users\danpo\Documents\MATLAB\ephysSuite\zf son mdy';
    db = containers.Map; % initialized as none, basically
end
```

Next, it's a simple matter of initializing your dbHandler!

```Matlab
dbh = dhHandler();
```
See the `@dbHandler` folder to learn more about how to use this class.
