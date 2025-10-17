#!/bin/env python
##
# @file This file is part of the ExaHyPE project.
# @author ExaHyPE Group (exahype@lists.lrz.de)
#
# @section LICENSE
#
# Copyright (c) 2016  http://exahype.eu
# All rights reserved.
#
# The project has received funding from the European Union's Horizon
# 2020 research and innovation programme under grant agreement
# No 671698. For copyrights and licensing, please consult the webpage.
#
# Released under the BSD 3 Open Source License.
# For the full license text, see LICENSE.txt
#
#
# @section DESCRIPTION
#
# Controller of the kernel generator
#
# @note
# requires python3


import os
import copy
import errno
import time

from .configuration import Configuration
from .argumentParser import ArgumentParser

from .models import *


class Controller:
    """Main Controller
    
    Read the input from the public API, validate them and generate a base 
    context for the models.
    
    Use generateCode() to run the models with the base context.
    
    Can generate gemms with generateGemms(outputFile, matmulconfig), will be done
    automatically when using generateCode().
    """
    
    def __init__(self, inputConfig = None):
        """Initialize the base config from the command line inputs"""
        if inputConfig == None:
            args = ArgumentParser.parseArgs()
        else:
            ArgumentParser.validateInputConfig(inputConfig)
            args = inputConfig
        
        self.commandLine = ArgumentParser.buildCommandLineFromConfig(args)
        
        # Generate the base config from the args input
        self.config = {
            "kernelType"            : args["kernelType"],
            "pathToOptKernel"       : args["pathToOptKernel"],
            "pathToOutputDirectory" : os.path.join(Configuration.pathToExaHyPERoot, args["pathToApplication"], args["pathToOptKernel"]),
            "solverName"            : args["solverName"],
            "codeNamespace"         : args["namespace"],
            "tempVarsOnStack"       : args["tempVarsOnStack"],
            "architecture"          : args["architecture"],
            "useLibxsmm"            : Configuration.matmulLib == "Libxsmm",
            "useEigen"              : Configuration.matmulLib == "Eigen",
            "pathToLibxsmmGemmGenerator"  : Configuration.pathToLibxsmmGemmGenerator,
            "runtimeDebug"          : Configuration.runtimeDebug, #for debug
            "prefetchInputs"        : Configuration.prefetching in ["Inputs", "All"],
            "prefetchOutputs"       : Configuration.prefetching in ["Outputs", "All"],
            "prefetchLevel"         : Configuration.prefetchLevel
        }
        
        if self.config["kernelType"] == "aderdg":
            self.config.update( {
                "numerics"              : args["numerics"],
                "nVar"                  : args["numberOfVariables"],
                "nPar"                  : args["numberOfParameters"],
                "nData"                 : args["numberOfVariables"] + args["numberOfParameters"],
                "nDof"                  : (args["order"])+1,
                "nDim"                  : args["dimension"],
                "useFlux"               : (args["useFlux"] or args["useFluxVect"]) or args["useViscousFlux"],
                "useFluxVect"           : args["useFluxVect"],
                "useViscousFlux"        : args["useViscousFlux"],
                "useNCP"                : (args["useNCP"] or args["useNCPVect"]),
                "useNCPVect"            : args["useNCPVect"],
                "useSource"             : (args["useSource"] or args["useSourceVect"] or args["useFusedSource"] or args["useFusedSourceVect"]),
                "useSourceVect"         : args["useSourceVect"],
                "useFusedSource"        : (args["useFusedSource"] or args["useFusedSourceVect"]),
                "useFusedSourceVect"    : args["useFusedSourceVect"],
                "nPointSources"         : args["usePointSources"],
                "usePointSources"       : args["usePointSources"] >= 0,
                "useMaterialParam"      : (args["useMaterialParam"] or args["useMaterialParamVect"]),
                "useMaterialParamVect"  : args["useMaterialParamVect"],
                "quadratureType"        : ("Gauss-Lobatto" if args["useGaussLobatto"] else "Gauss-Legendre"),
                "useCERKGuess"          : args["useCERKGuess"],
                "useSplitCK"            : args["useSplitCK"],
                "useVectPDE"            : args["useVectPDE"],
                "useAoSoA2"             : args["useAoSoA2"],
                "predictorRecompute"    : args["predictorRecompute"],
                "initialGuess"          : "mixedPicard" #TODO JMG put as proper toolkit arg
                #"initialGuess"          : "default" #TODO JMG put as proper toolkit arg
            })
            self.config["useSourceOrNCP"] = self.config["useSource"] or self.config["useNCP"]
        elif self.config["kernelType"] == "limiter":
            self.config.update( {
                "nVar"                  : args["numberOfVariables"],
                "nPar"                  : args["numberOfParameters"],
                "nData"                 : args["numberOfVariables"] + args["numberOfParameters"],
                "nDof"                  : (args["order"])+1,
                "nDofLim"               : args["limPatchSize"] if args["limPatchSize"] >=0 else 2*args["order"]+1,
                "nDim"                  : args["dimension"],
                "nObs"                  : args["numberOfObservable"],
                "ghostLayerWidth"       : args["ghostLayerWidth"],
                "quadratureType"        : ("Gauss-Lobatto" if args["useGaussLobatto"] else "Gauss-Legendre")
            })
        elif self.config["kernelType"] == "fv":
            self.config.update( {
                "nVar"                  : args["numberOfVariables"],
                "nPar"                  : args["numberOfParameters"],
                "nData"                 : args["numberOfVariables"] + args["numberOfParameters"],
                "nDof"                  : args["patchSize"],
                "nDim"                  : args["dimension"],
                "useFlux"               : args["useFlux"] or args["useViscousFlux"],
                "useViscousFlux"        : args["useViscousFlux"],
                "useNCP"                : args["useNCP"],
                "useSource"             : args["useSource"] or args["useFusedSource"],
                "useFusedSource"        : args["useFusedSource"],
                "useRobustDiagLim"      : args["useRobustDiagLim"],
                "nPointSources"         : args["usePointSources"],
                "usePointSources"       : args["usePointSources"] >= 0,
                "useMaterialParam"      : args["useMaterialParam"],
                "finiteVolumesType"     : args["finiteVolumesType"],
                "ghostLayerWidth"       : 2, #hard coded musclhancock value
                "slopeLimiter"          : args["slopeLimiter"]
            })
            self.config["useSourceOrNCP"] = self.config["useSource"] or self.config["useNCP"]
            
        self.validateConfig(Configuration.simdWidth.keys())
        self.config["vectSize"] = Configuration.simdWidth[self.config["architecture"]] #only initialize once architecture has been validated
        self.config["cachelineSize"] = Configuration.cachelineSize[self.config["architecture"]] #only initialize once architecture has been validated
        self.baseContext = self.generateBaseContext() # default context build from config
        self.matmulList = [] #list to store the tupple (fileName, matmulConfig) of all requested Matmul (used for gemmsGeneratorModel)
        self.gemmList = [] #list to store the name of all generated gemms (used for gemmsCPPModel)

    def validateConfig(self, validArchitectures):
        """Ensure the configuration fit some constraint, raise errors if not"""
        if not (self.config["architecture"] in validArchitectures):
           raise ValueError("Architecture not recognized. Available architecture: "+str(validArchitectures))
        if self.config["kernelType"] == "aderdg":
            if not (self.config["numerics"] == "linear" or self.config["numerics"] == "nonlinear"):
                raise ValueError("numerics has to be linear or nonlinear")
            if self.config["nVar"] < 0:
               raise ValueError("Number of variables must be >=0 ")
            if self.config["nPar"] < 0:
               raise ValueError("Number of parameters must be >= 0")
            if self.config["nDim"] < 2 or self.config["nDim"] > 3:
               raise ValueError("Number of dimensions must be 2 or 3")
            if self.config["nDof"] < 1 or self.config["nDof"] > 15+1: #nDof = order+1
               raise ValueError("Order has to be between 0 and 15 (inclusive)")
        if self.config["kernelType"] == "fv":
            if self.config["finiteVolumesType"] != "musclhancock":
               raise ValueError("Only musclhancock scheme is supported")


    def printConfig(self):
        print(self.config)


    def generateBaseContext(self):
        """Generate a base context for the models from the config (use hard copy)"""
        context = copy.copy(self.config)
        context["nVarPad"]  = self.getSizeWithPadding(context["nVar"])
        context["nParPad"]  = self.getSizeWithPadding(context["nPar"])
        context["nDataPad"] = self.getSizeWithPadding(context["nData"])
        context["nDofPad"]  = self.getSizeWithPadding(context["nDof"])
        context["nDof3D"]   = 1 if context["nDim"] == 2 else context["nDof"]
        context["nDof2Pad"] = self.getSizeWithPadding(context["nDof"]*context["nDof"])
        context["solverHeader"]      = context["solverName"].split("::")[1] + ".h"
        context["codeNamespaceList"] = context["codeNamespace"].split("::")
        context["guardNamespace"]    = "_".join(context["codeNamespaceList"]).upper()
        if self.config["kernelType"] == "limiter":
            context["nDofLimPad"] = self.getSizeWithPadding(context["nDofLim"])
            context["nDofLim3D"] = 1 if context["nDim"] == 2 else context["nDofLim"]
            context["ghostLayerWidth3D"] = 0 if context["nDim"] == 2 else context["ghostLayerWidth"]
        elif self.config["kernelType"] == "aderdg":
            context["isLinear"] = context["numerics"] == "linear"
            context["useVectPDEs"] = context["useFluxVect"] or True #TODO JMG add other vect
        elif self.config["kernelType"] == "fv":
            context["ghostLayerWidth3D"] = 0 if context["nDim"] == 2 else context["ghostLayerWidth"]
            context["nDofG"] = context["ghostLayerWidth"]*2 + context["nDof"]
            context["nDofG3D"] = 1 if context["nDim"] == 2 else context["nDofG"]
            
        context["function_utils"] = {"getPadSize": self.getPadSize} # quadratureModel needs access to the controller function to calculate padding size
        return context

    def getSizeWithPadding(self, sizeWithoutPadding):
        """Return the size of the input with the architecture specific padding added"""
        return self.config["vectSize"] * int((sizeWithoutPadding+(self.config["vectSize"]-1))/self.config["vectSize"])


    def getPadSize(self, sizeWithoutPadding):
        """Return the size of padding required for its input"""
        return self.getSizeWithPadding(sizeWithoutPadding) - sizeWithoutPadding


    def generateCode(self):
        """Main method: call the models to generate the code"""
        
        # create directory for output files if not existing
        try:
            os.makedirs(self.config["pathToOutputDirectory"])
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise
        
        # remove all .cpp, .cpph, .c and .h files (we are in append mode!) as well as previous symlink (see symlinkBLASlib())
        for fileName in os.listdir(self.config["pathToOutputDirectory"]):
            _ , ext = os.path.splitext(fileName)
            if ext in [".cpp", ".cpph", ".c", ".h"] or fileName in ["Eigen"]:
                os.remove(self.config["pathToOutputDirectory"] + "/" + fileName)
        
        # Symlink the BLAS library if needed
        self.symlinkBLASlib()
        
        # run the models new files
        
        self.runModel(    "kernelsHeader",            kernelsHeaderModel.KernelsHeaderModel(self.baseContext))
        self.runModel(    "configurationParameters",  configurationParametersModel.ConfigurationParametersModel(self.baseContext))
        
        if self.config["kernelType"] in ["aderdg", "limiter"]:
            self.runModel("quadrature",               quadratureModel.QuadratureModel(self.baseContext))
        
        if self.config["kernelType"] == "aderdg":
            self.runModel("converter",                converterModel.ConverterModel(self.baseContext))
            self.runModel("amrRoutines",              amrRoutinesModel.AMRRoutinesModel(self.baseContext))
            self.runModel("deltaDistribution",        deltaDistributionModel.DeltaDistributionModel(self.baseContext))
            self.runModel("faceIntegral",             faceIntegralModel.FaceIntegralModel(self.baseContext))
            self.runModel("fusedSTPVI",               fusedSpaceTimePredictorVolumeIntegralModel.FusedSpaceTimePredictorVolumeIntegralModel(self.baseContext))
            self.runModel("matrixUtils",              matrixUtilsModel.MatrixUtilsModel(self.baseContext))
            self.runModel("dgMatrix",                 dgMatrixModel.DGMatrixModel(self.baseContext))
            self.runModel("solutionUpdate",           solutionUpdateModel.SolutionUpdateModel(self.baseContext))
            self.runModel("surfaceIntegral",          surfaceIntegralModel.SurfaceIntegralModel(self.baseContext))
        
        if self.config["kernelType"] == "limiter":
            self.runModel("limiter",                  limiterModel.LimiterModel(self.baseContext, self))
        
        if self.config["kernelType"] == "fv":
            self.runModel("ghostLayerFilling",        fvGhostLayerFillingModel.FVGhostLayerFillingModel(self.baseContext))
            self.runModel("ghostLayerFillingAtBoundary", fvGhostLayerFillingAtBoundaryModel.FVGhostLayerFillingAtBoundaryModel(self.baseContext))
            self.runModel("boundaryLayerExtraction",  fvBoundaryLayerExtractionModel.FVBoundaryLayerExtractionModel(self.baseContext))
            self.runModel("solutionUpdate",           fvSolutionUpdateModel.FVSolutionUpdateModel(self.baseContext))
        
        if self.config["kernelType"] in ["aderdg", "fv"]:
            self.runModel("boundaryConditions",       boundaryConditionsModel.BoundaryConditionsModel(self.baseContext))
            self.runModel("stableTimeStepSize",       stableTimeStepSizeModel.StableTimeStepSizeModel(self.baseContext))
            self.runModel("adjustSolution",           adjustSolutionModel.AdjustSolutionModel(self.baseContext))
            self.runModel("riemannSolver",            riemannModel.RiemannModel(self.baseContext))
        
        ## must be run only after all gemm's configurations have been generated
        gemmsGeneratorContext = copy.copy(self.baseContext)
        gemmsGeneratorContext["matmulList"] = self.matmulList
        self.runModel(    "gemmsGenerator",           gemmsGeneratorModel.GemmsGeneratorModel(gemmsGeneratorContext))
        gemmsContext = copy.copy(self.baseContext)
        gemmsContext["gemmList"] = self.gemmList
        self.runModel(    "gemmsCPP",                 gemmsCPPModel.GemmsCPPModel(gemmsContext))


    def runModel(self, name, model):
        """Run the given model and if debug then print runtime"""
        start = time.perf_counter()
        out = model.generateCode()
        if out is not None and "matmulList" in out and "fileName" in out:
            self.matmulList += map(lambda m: (out["fileName"], m), out["matmulList"])
        if out is not None and "gemmList" in out:
            self.gemmList += out["gemmList"]
        if self.config["runtimeDebug"]:
            t = time.perf_counter() - start
            print(name+": "+str(t))

            
    def symlinkBLASlib(self):
        if self.config["useEigen"]:
            os.symlink(os.path.join(Configuration.pathToEigen, "Eigen"), os.path.join(self.config["pathToOutputDirectory"], "Eigen"))
