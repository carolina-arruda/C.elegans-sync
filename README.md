# C.elegans-sync
C. elegans neural network - synchronization and modularity

This repository shares all information used on this paper https://www.sciencedirect.com/science/article/pii/S0378437119311938.

In this work we study the modular structure of the C. elegans neural network using the electrical junction connections. The system is composed of 248 neurons (nodes) and 511 gap junctions (links). The adjacency matrix is symmetric and weighted. The complete Connectivity Data can be find at https://www.wormatlas.org/neuronalwiring.html#Connectivitydata (accessed March 2020).

Files:

1) Neuron1_Neuron2_Weight.csv: contains the connections between "Neuron 1" and "Neuron 2" and its respective number of gap junctions (column Weight). 

2) EJ248adj.csv: is the adjacency matrix.

3) Data-EJ248.csv: the columns specify the name of Neuron, the Index label, physical Localization	(head, mid-body, tail), the Topology classification into 3 modules, the functional Class	(interneuron, sensory neurons and moroneurons), the type of Ganglion (A, B, C, D, E, F, G, H, J, K) and the Topology classification into 5 modules and 10 modules. It also contains the number of neurons in each group used.

The codes used to run the partially forced Kuramoto model and the correlation matrices are kuramoto_celegans.f90 and correlation_celegans.f90. Type gfortran program.f90 to compile and ./a.out to run.
