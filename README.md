# Effects of transcranial alternating current stimulation on spiking activity in computational models of single neocortical neurons

This folder contains the codes used in the paper Tran et al., 2020. Effects of transcranial alternating current stimulation on spiking activity in computational models of single neocortical neurons. In this paper, we investigated the effects of tACS amplitudes on computational models of single neocortical neurons.

## Information about the data
Due to the high size of data, not all of them will be available in that depository. Some are given and will be used as an example. The different neurons model used are those used in aberra2018biophysically.



## Data generation
The data were generated using NEURON environement 7.8.1  You can download this software here https://neuron.yale.edu/neuron/. A supercomputer was used, however, depending on the length of simulation, it is possible to use a regular desktop or laptop. For sub(supra)threshold investigation, the code is quite similar. For the suprathreshold activity investigation, a synapse was modeled. If you are investigating the subthreshold activity, a silent neuron is used (no synapse). In that case, please comment line 92 from the init.hoc file.

## Data analysis
To analyse the data, a custom matlab code was used.



## Reference
Please use the following reference: Tran et al., 2020. Effects of transcranial alternating current stimulation on spiking activity in computational models of single neocortical neurons

## Correspondence
htran - at- umn.edu (Harry Tran)

aopitz -at- umn.edu (Alex Opitz)
