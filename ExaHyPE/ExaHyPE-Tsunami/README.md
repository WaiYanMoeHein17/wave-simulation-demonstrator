# This is the ExaHyPE project #

This repository mainly exists to run the Tohoku tsunami scenario via https://hub.docker.com/r/linusseelinger/exahype.

For any other applications or uses please use the main ExaHyPE repository at: https://gitlab.lrz.de/exahype/ExaHyPE-Engine

### It contains the following directories:
ApplicationExamples: Various refinements of the Tohoku-tsunami scenario
CodeGenerator: The optimized kernels generators
Peano: Mounting point for the Peano installation
Submodules: Installation directory for the git submodules (ExaHyPE's dependencies)
Toolkit: The ExaHyPE development toolkit

## Mini installation guide
No installation necessary, can be run directly through docker via https://hub.docker.com/r/linusseelinger/exahype.

### Commit guidelines
Please, don't commit the following:

Binary files (*.o, executables, ... ) excluding those necessary for the documentation
Output files (*.vtk, logs, ... )

This project has received funding from the European Union’s Horizon 2020 research and innovation programme under grant agreement No 671698.
# ExaHyPE-Sandbox
