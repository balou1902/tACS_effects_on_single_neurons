# Effects of transcranial alternating current stimulation on spiking activity in computational models of single neocortical neurons

This folder contains the codes used in the paper Tran et al., 2020. Effects of transcranial alternating current stimulation on spiking activity in computational models of single neocortical neurons. In this paper, we investigated the effects of tACS amplitudes on computational models of single neocortical neurons.

## Information about the data
A set of 25 neocortical cells was used in this investigation and can be downloaded from ModelDB (accession number: 066023 or directly click on this link https://senselab.med.yale.edu/modeldb/ShowModel?model=241165#tabs-1).
Further information about the cell dynamics can be found in [aberra2018biophysically].


## Data generation
The data were generated using NEURON environement 7.8.1  You can download this software here https://neuron.yale.edu/neuron/. A supercomputer was used, however, depending on the length of simulation, it is possible to use a regular desktop or laptop. For sub(supra)threshold investigation, the code is quite similar. For the suprathreshold activity investigation, a synapse was modeled. If you are investigating the subthreshold activity, a silent neuron is used (no synapse). In that case, please comment line 92 from the init.hoc file.

## Data analysis
To analyse the data, a custom matlab code was used.

### Subthreshold activity analysis
Run the code 'subthreshold_analysis.m'. This code analyzes the somatic membrane fluctuations for silent neurons and computes the polarization length for each cell. This code generates each of the subplots of the figure 2 of the article. Please notice that a path change for the matlab working directory is required.

### Suprathreshold activity analysis



## Reference
Please use the following reference: Tran et al., 2020. Effects of transcranial alternating current stimulation on spiking activity in computational models of single neocortical neurons

## Correspondence
htran - at- umn.edu (Harry Tran)

aopitz -at- umn.edu (Alex Opitz)
