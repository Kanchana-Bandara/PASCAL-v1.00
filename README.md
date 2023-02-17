## Pan-Arctic Behavioural and Life-history Simulator, PASCAL (v1.00)
#### Generalized high-latitude copepod simulation model (baseline)
This repository contains the first version of PASCAL, initiated in 2016 & completed in 2018. PASCAL v.1.00 contains a general (i.e., non-species-specific), high-resolution behavioural and life-history simulation model for a generalized copepod resembling a synthetic _C. finmarchicus_ - _C. glacialis_ hybrid.

The outputs from this model have already been analyzed and published once at:

_Bandara, K., Varpe, Ø., Ji, R., & Eiane, K. (2018). A high-resolution modeling study on diel and seasonal vertical migrations of high-latitude copepods. Ecological Modelling, 368, 357-376._

##### Contents
This repository contains two main files: (i) the PASCAL v1.00 hull and (ii) PASCAL v1.00 simulation drive. The hull file is the main file from which the simulation function is called. The additional functions file contains three (3) small functions to be loaded besides the two main files.

##### Warnings
The environment files (2D arrays) are needed to execute the model. These files are not included here. Check the corresponding Zenodo repository for example files.
The HPC framework has been redarcted due to an ongoing publication. We will add this as soon as we can in a separate repository (this README file will be updated when the HPC framework is made Open Access)

###### Runtime
The model is set-up to run for 100 generations (calendar years). With the HPC framework built-in, it takes about 18-21 hours to complete execution in a CORE i9 7920X 24-Thread, 2.9-3.5 GHz gaming rig. If you have the PushoverR app, the model will send an update to your phone once the simulation has been completed. For single threaded operations (no HPC), please dial down the no. of individuals in the model. However, the impact of the simulation population size on the emergent ecological properties has not been assessed. 

##### Lisence
CC BY 4.0
