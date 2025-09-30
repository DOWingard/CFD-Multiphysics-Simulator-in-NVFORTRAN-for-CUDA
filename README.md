# CFD MultiPhysics Finite Element Simulator in NVFORTRAN for parallelized CPU or GPU
Currently just simulates water with a pressure gradient and suffers to be numerically stable

## Current Simulation
The manifold is a 3d mesh of cell objects which can be assigned a conserved quantity vector U(:) and uses a subroutine to set it up from water. Will me modularized later impliment a standard process of adding a different simiulation kernel (ie plasma). The manifold gets updated by calculating a tangent bundle, where the Rusanov flux across every face is calculated and stored in tangent bundle mesh which can be used to write the t+1 state of the corresponding cell of the main manifold. The comparison of t,t+1 values is used for smoothing to help numerical stability (doesn't do a lot lol).

Simulation runs and logs time.

**Current Setup**
Runs a 40x40x40 Finite Element model over .2 seconds (100 timesteps of 0.002s ) 
* gets pretty unstable quickly after this
* computes on parallelized CPU threads (12 threads with 16GB RAM)

## Working on the Sim
I have some .vscode commands to launch the container, assuming you have docker desktop (at least CLI) set up and the CUDA compiler files set up in /opt/ in sudo directory. 

**Launch Container**

run task Start Docker Container

**Set Up Environment**

run task Setup Environment

**Compiling Fortran**
There is a makefile system set up so simply 
```bash
cd kernel
```
then you can do 
```bash
make
```
to compil and 
```bash
make clean
```
to clear the build files. 
```bash
empty
```
clears all the data/ folder files. I also have gdb set up so you can run 
```bash
make debug
```
and use gdb to debug it. 



## Scripts 
I have some QOL scripts

__in ROOT__

* gitmain   -> git push to main branch
* gitster   -> "" master
* gitbranch -> "" $USER_INPUT_BRANCH

__in ROOT/kernel__

* empty -> cleans data/ folder


## Scalability
The model gets allocated on host and, while not modularized at the moment, one can add a different kernel, cell computations, and tangent bundle tensor to change the simulation 

**Important TODOs**

* Allow all the properties written to the conserved quanitities vector to be accessed dynamically so feedback can be incorporated for multiphysics simluation
* MAKE MORE NUMERICALLY STABLE
    * Use chebyshev routine to calculate differentials using recursion where possible
* GPU Methods
* Pipeline the build process that wraps with python to plot outputs


