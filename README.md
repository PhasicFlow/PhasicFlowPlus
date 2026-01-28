<div align ="center">
<img src="doc/phasicFlowPlus_Logo_github.png" style="width: 400px;">
</div>

.

**PhasicFlowPlus** is a software package for simulating fluid-particle flows. It is a combination of computational fluid dynamics (CFD) and discrete element method (DEM). The fluid is assumed to be continuum phase and the particles as discrete bodies.

Here, DEM calculations are handled using features of [PhasicFlow](https://github.com/PhasicFlow/phasicFlow). It is a parallel DEM package that can be run on multi-core CPUs or GPUs. The equations for the fluid phase are discretized based on finite volume (FV) and solved using OpenFOAM. OpenFOAM is parallelized based on MPI for running on on multicore CPUs. The fluid-particle coupling (PhasicFlowPlus) uses both parallelization methods, shared-memory and MPI, to leverage the maximum computational resoureces. 

Based on the above configuration, PhasicFlowPlus can use the computational resources of a multi-core CPU or use the computational resource of both CPU and GPU. 

## What is under development?
The following parts are being developed at the moment:
* Coarse-graining for CFD-DEM
* Resolved solver for CFD-DEM
* Modifying some parts for better functionality and performance  

## How to build PhasicFlowPlus
### Version-1.0
Version-1.0 is compatible with PhasicFlow-v-1.0 and OpenFOAM-v24. You need to [install phasicFlow-v-1.0](https://github.com/PhasicFlow/phasicFlow/wiki/How-to-build-PhasicFlow%E2%80%90v%E2%80%901.0) and then [install OpenFOAM-v24](https://www.cemf.ir/installing-openfoam-v2412-on-ubuntu-and-windows/) on your computer. After that, follow the following steps to install PhasicFlowPlus.

+ **Step 1:** copy the source code into `~/PhasicFlow` folder:

    ```bash
    cd ~/PhasicFlow
    wget https://github.com/PhasicFlow/PhasicFlowPlus/archive/refs/heads/main.zip
    unzip -q main.zip
    mv PhasicFlowPlus-main/ PhasicFlowPlus
    # Note:
    #   Instead of using wget, you could directly clone it in director ~/PhasicFlow
    #   use the command: git clone https://github.com/PhasicFlow/PhasicFlowPlus.git  
    ```
    
+ **Step 2:** build the software:
    ```bash
    # Note:
    #   make sure that the phasicFlow is correctly installed with checkPhasicFlow
    #   First activate OpenFOAM-v2412 (choose the command that matches your installation)
    #   1- pre-built: source /usr/lib/openfoam/openfoam2412/etc/bashrc
    #   2- source pack : source $HOME/OpenFOAM-v2412/OpenFOAM-v2412/etc/bashrc
    cd ~/PhasicFlow/PhasicFlowPlus
    ./Allwmake
    ```

### Version-0.1 
First, you need to [install PhasicFlow-v0.1](https://github.com/PhasicFlow/phasicFlow/wiki/How-to-Build-PhasicFlow) and OpenFoam-v9 (For now, it is only tested with [OpenFOAM-v9](https://openfoam.org/download/9-source/)) on your computer. After that, copy [PhasicFlowPlus-v0.1](https://github.com/PhasicFlow/PhasicFlowPlus/releases/tag/v-0.1) on your computer. The `PhasicFlowPlus` folder should be located beside `phasicFlow` folder on your computer (in ~/PhasicFlow/ folder). Navigate to the root directory of the code and enter the following commands to install the code.

```bash
# Note:
#   make sure that the phasicFlow is correctly installed with checkPhasicFlow
#   First activate OpenFOAM-v9 (choose the command that matches your installation)
#   1- pre-built: source /opt/openfoam9/etc/bashrc
#   2- source pack : source $HOME/OpenFOAM/OpenFOAM-9/etc/bashrc
cd ~/PhasicFlow/PhasicFlowPlus/ 
./Allwmake
```

## Core Developers and Contributors

PhasicFlowPlus is developed and maintained by a dedicated team of researchers and engineers specializing in computational fluid dynamics (CFD) and discrete element method (DEM). The core developers include:

* Hamidreza Norouzi
* Bahram Haddadi

Contributions from the community are highly encouraged. Contributors can submit bug fixes, new features, or documentation improvements through pull requests on the GitHub repository. In the current version, we also had contributors from (the list will be updated):

* Alireza Hosseini
* Nima Joghataei
