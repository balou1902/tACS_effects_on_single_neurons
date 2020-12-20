
# NEURON Model generation

This folder contains all the required .hoc files to generate the somatic membrane voltage for the cells used in the investigation.

## Files description

### init.hoc

### mySettings.hoc

In this file, the user enters the main parameters of the stimulation (stimulation length and frequency) as well as external paramaters: the index of the cell used for the simulation, the orientation of the electric field and the use of a synaptic input.

### mySynapse.hoc

This files generates the synaptic input and the user can modify parameters such as the synaptic location, the type of synaptic modeling or the connection parameters.

### exportLocs.hoc

### stimwaveform.hoc
This file generates the waveform of the TACS using the different inputs (stimulation frequency, amplitudes, stimulation length).



