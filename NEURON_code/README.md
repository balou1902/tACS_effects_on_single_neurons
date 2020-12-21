
# NEURON Model generation

This folder contains all the required .hoc files to generate the somatic membrane voltage for the cells used in the investigation.

## Files description

### init.hoc

This file is the main file of the simulation. All the .hoc files required to run a simulation are loaded. The user can use a simple value for the tACS amplitude or a vector of values. At the end of the simulation, a folder will be generated in which results will be saved: locations of all segments, somatic membrane voltage, the tACS waveform and the time vector. This is the file you need to run to launch a simulation.
Please notice that you need to modify the path of your working directory.

### mySettings.hoc

In this file, the user enters the main parameters of the stimulation (stimulation length and frequency) as well as external paramaters: the type of cell used for the simulation, the orientation of the electric field and the use of a synaptic input.

### mySynapse.hoc

This files generates the synaptic input and the user can modify parameters such as the synaptic location, the type of synaptic modeling or the connection parameters.

### exportLocs.hoc

This file saves all the coordinates of each segment of a neuron morphology.

### stimwaveform.hoc
This file generates the waveform of the TACS using the different inputs (stimulation frequency, amplitudes, stimulation length).



