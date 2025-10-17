Prerequisites
=============

Python
------

Python 3.3 or more is required.

Dependencies
------------

By default get its dependencies from ExaHyPE's submodules

You can install locally all the dependencies using the provided bash script:

`./KernelGenerator/importJinjaAndUseSource.sh`

You then need to adapt `./KernelGenerator/kernelgenerator/configuration.py` with 
the correct paths to the dependencies


### Jinja2

The KernelGenerator uses the template engine Jinja2 (http://jinja.pocoo.org/)

If it is not provided by your python installation you can instead use the source 
of jinja2 directly either with the bash script or by following these steps

1) Clone the source from the git repository to your desired install directory <my-path>: 
`git clone https://github.com/pallets/jinja.git <my-path>/jinja`

2) In `ExaHyPE-Engine/KernelGenerator/kernelgenerator/configuration.py`, 
put the correct path to the jinja directory

3) Do the same for Markupsafe (https://pypi.org/project/MarkupSafe/), jinja's
dependency


### LIBXSMM

The CodeGenerator uses LIBXSMM to perform efficient matrux-matrix operations

You can either install it with the bash script or by following these steps

1) Clone the source from the git repository to your desired install directory <my-path>: 
`git clone -b release --single-branch https://github.com/hfp/libxsmm.git <my-path>/libxsmm`

2) Compile the source code with 
`make realclean && make generator`

3) In `ExaHyPE-Engine/KernelGenerator/kernelgenerator/configuration.py`, 
put the correct path to the libxsmm\_gemm\_generator


Paths
-----

Every path is relative to the root of the project (inside the ExaHyPE-Engine 
directory).

The KernelGenerator relies on the paths provided in `configuration.py`:

The generated code will be put accordingly to the `pathToOptKernel` argument 
starting from the internal ExaHyPe, by default in a `kernel` subdirectory of your 
application.


KernelGenerator
===============

The KernelGenerator can be used by other python program using it's python3 API or
be called directly using the command line.

It is automatically called by ExaHyPE's Toolkit through the python3 API if 
optimised kernels are specified.
A verbose Toolkit displays the command lines equivalent to its API calls.

Using the command line interface, to access the global help use: 
`python3 KernelGenerator/kernelgenerator -h`

To access the kernel specific help use: 
* ADERDG:        `python3 KernelGenerator/kernelgenerator aderdg -h`
* Finite Volume: `python3 KernelGenerator/kernelgenerator fv -h`
* Limiter:       `python3 KernelGenerator/kernelgenerator limiter -h`


ADERDG -  usage and arguments
-----------------------------

positional arguments:
*  kernelType            aderdg
*  pathToApplication     path to the application as given by the ExaHyPE specification file (application directory as root)
*  pathToOptKernel       desired relative path to the generated code (application directory as root)
*  namespace             desired namespace for the generated code
*  solverName            name of the user-solver
*  numberOfVariables     the number of quantities
*  numberOfParameters    the number of parameters (fixed quantities)
*  order                 the order of the approximation polynomial
*  dimension             the number of spatial dimensions in the simulation (2 or 3)
*  numerics              linear or nonlinear
*  architecture          the microarchitecture of the target device

optional arguments:
*  -h, --help            show this help message and exit
*  --useFlux             enable flux
*  --useFluxVect         enable vectorized flux (include useFlux)
*  --useViscousFlux      enable viscous flux
*  --useNCP              enable non conservative product
*  --useNCPVect          enable vectorized non conservative product (include useNCP)
*  --useSource           enable source terms
*  --useSourceVect       enable vectorized source terms (include useSource)
*  --useFusedSource      enable fused source terms (include useSource)
*  --useFusedSourceVect  enable vectorized fused source terms (include useFusedSource and useSourceVect)
*  --useMaterialParam    enable material parameters
*  --useMaterialParamVect enable vectorized material parameters
*  --usePointSources numberOfPointSources enable numberOfPointSources point sources
*  --useCERKGuess        use CERK for SpaceTimePredictor inital guess (nonlinear only)
*  --useSplitCK          use split Cauchy–Kowalevski formulation (linear only)
*  --useGaussLobatto     use Gauss Lobatto Quadrature instead of Gauss Legendre
*  --useVectPDE          use vectorized PDE terms (applies when present to: Flux, NCP, Source, FusedSource and MaterialParam)
*  --tempVarsOnStack     put the big scratch arrays on the stack instead of the
                         heap (you can use ulimit -s to increase the stack size)

Example: ``python3 KernelGenerator/kernelgenerator aderdg nonlinear_Euler_Sod_3d kernels/Euler_EulerSolver/aderdg Euler::EulerSolver_ADERDG_kernels::aderdg Euler::EulerSolver 5 0 9 2 nonlinear hsw --useFlux --tempVarsOnStack``

Finite Volume - usage and arguments
-----------------------------------

positional arguments:
*  kernelType            fv
*  pathToApplication     path to the application as given by the ExaHyPE specification file (application directory as root)
*  pathToOptKernel       desired relative path to the generated code (application directory as root)
*  namespace             desired namespace for the generated code
*  solverName            name of the user-solver
*  numberOfVariables     the number of quantities
*  numberOfParameters    the number of parameters (fixed quantities)
*  patchSize             the size of a patch
*  dimension             the number of spatial dimensions in the simulation (2 or 3)
*  finiteVolumesType     linear or nonlinear
*  architecture          the microarchitecture of the target device
*  {minmod,koren,superbee,vanalbada,mclim} slope limiter function for the scheme

optional arguments:
*  -h, --help            show this help message and exit
*  --useFlux             enable flux
*  --useViscousFlux      enable viscous flux
*  --useNCP              enable non conservative product
*  --useSource           enable source terms
*  --useFusedSource      enable fused source terms (include useSource)
*  --useMaterialParam    enable material parameters
*  --useRobustDiagLim    enable robust diagonal limiting in musclhancock scheme
*  --usePointSources numberOfPointSources enable numberOfPointSources point sources
*  --tempVarsOnStack     put the big scratch arrays on the stack instead of the heap (you can use ulimit -s to increase the stack size)


Example: ``python3 KernelGenerator/kernelgenerator fv nonlinear_Euler_Sod_3d kernels/Euler_EulerSolver/fv Euler::EulerSolver_FV_kernels::fv Euler::EulerSolver 5 0 19 2 musclhancock hsw --useFlux``


Limiter - usage and arguments
-----------------------------

positional arguments:
*  kernelType            limiter
*  pathToApplication     path to the application as given by the ExaHyPE specification file (application directory as root)
*  pathToOptKernel       desired relative path to the generated code (application directory as root)
*  namespace             desired namespace for the generated code
*  solverName            name of the user-solver
*  numberOfVariables     the number of quantities
*  numberOfParameters    the number of parameters (fixed quantities)
*  order                 the order of the approximation polynomial
*  dimension             the number of spatial dimensions in the simulation (2 or 3)
*  architecture          the microarchitecture of the target device
*  numberOfObservable    number of observable
*  limPatchSize          the size of the limiter patches per coordinate direction.

optional arguments:
*  -h, --help            show this help message and exit
*  --ghostLayerWidth width use limiter with the given ghostLayerWidth, default = 0
*  --useGaussLobatto     use Gauss Lobatto Quadrature instead of Gauss Legendre
*  --tempVarsOnStack     put the big scratch arrays on the stack instead of the heap (you can use ulimit -s to increase the stack size)


Example: ``python3 KernelGenerator/kernelgenerator limiter nonlinear_Euler_Sod_3d kernels/Euler_EulerSolver/limiter Euler::EulerSolver_kernels::limiter Euler::EulerSolver 5 0 9 2 hsw 2 19 --ghostLayerWidth 2``


Data format and padding
-----------------------

The KernelGenerator may use padding when producing architecture specific code, 
it may also change the index order

Conversion between the generic and optimised format can be done using the 
generated converter

Using the C index order convention with index in italic being absent in dim 2

Advanced options like ``--useSplitCK`` or ``--useVectPDE`` may use other data layout, see the 
respective templates


| Array | Generic | Optimised | Note |
| ----- | ------- | --------- | ---- | 
| luh | _nDof_, nDof, nDof, nData | _nDof_, nDof, nDof, nData | unchanged for compatibility purpose|
| lQhbnd | 2*nDim, _nDof_, nDof, **nData** | 2*nDim, _nDof_, nDof, **nDataPad** | 2*nDim face on the square/cube |
| lFhbnd |  2*nDim, _nDof_, nDof, **nVar** | 2*nDim, _nDof_, nDof, **nVarPad** | 2*nDim face on the square/cube |
| lQhi |  _nDof_, nDof, nDof, **nData** | _nDof_, nDof, nDof, **nDataPad** | |
| lFhi |  nDim+1, _nDof_, nDof, nDof, **nVar** | nDim+1, _nDof_, nDof, nDof, **nVarPad** | lFhi has nDim+1 blocks == each directions + source |
| LQi **NL** |  nDof, _nDof_, nDof, nDof, **nData** | _nDof_, nDof, nDof, nDof, **nDataPad** | nonlinear case |
| LQi **Lin** |  nDof+1, _nDof_, nDof, nDof, **nData** | nDof+1, _nDof_, nDof, nDof, **nDataPad** | linear case |
| LFi **NL** |  nDof, _nDof_, nDof, nDof, nDim+1, **nVar** | nDim+1, nDof, _nDof_, nDof, nDof, **nVarPad** | nonlinear case, +1 for source. Move dimension coordinate to mimick lFhi blocks |
| LFi **Lin** |  2*nDim+1, nDof, _nDof_, nDof, nDof, **nVar** | 2*nDim+1, nDof, _nDof_, nDof, nDof, **nVarPad** | linear case, +1 for source. |
| rhs |  _nDof_, nDof, nDof, nDof, **nData** | _nDof_, nDof, nDof, nDof, **nDataPad** | nonlinear case |
| gradQ |  _nDof_, nDof, nDof, nDof, nDim, **nVar** | _nDof_, nDof, nDof, nDof, nDim, **nVarPad** | Not used if no NCP |
