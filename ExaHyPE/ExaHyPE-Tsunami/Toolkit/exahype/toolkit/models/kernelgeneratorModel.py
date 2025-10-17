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
# Generates the application-specific Makefile
#


import sys

# add path to dependencies
from ..configuration import Configuration
sys.path.insert(1, Configuration.pathToKernelgenerator)
import kernelgenerator


class KernelgeneratorModel:


    def generateCode(self, solverContext):
        #translation from solverContext to kernelgeneratorContext
        if solverContext["kernelType"] == "aderdg":
            kernelgeneratorContext = {
                # Mandatory parameters
                "kernelType"         : "aderdg",
                "pathToApplication"  : solverContext["outputPath"],
                "pathToOptKernel"    : solverContext["optKernelPath"],
                "namespace"          : solverContext["optNamespace"],
                "solverName"         : solverContext["project"] + "::" + solverContext["solver"],
                "numberOfVariables"  : solverContext["numberOfVariables"],
                "numberOfParameters" : solverContext["numberOfMaterialParameters"],
                "order"              : solverContext["order"],
                "dimension"          : solverContext["dimensions"],
                "numerics"           : "linear" if solverContext["isLinear"] else "nonlinear",
                "architecture"       : solverContext["architecture"],
                # Optional bool parameters (may set redundant flags and default false flag)
                "useFlux"            : solverContext["useFlux"],
                "useFluxVect"        : solverContext["useFluxVect"], #TODO JMG: legacy, use usePDEVect instead
                "useViscousFlux"     : solverContext["useViscousFlux"],
                "useNCP"             : solverContext["useNCP"],
                "useNCPVect"         : solverContext["useNCPVect"], #TODO JMG: legacy, use usePDEVect instead
                "useSource"          : solverContext["useSource"],
                "useSourceVect"      : solverContext["useSourceVect"],  #TODO JMG: legacy, use usePDEVect instead
                "useFusedSource"     : solverContext["useFusedSource"],
                "useFusedSourceVect" : solverContext["useFusedSourceVect"],  #TODO JMG: legacy, use usePDEVect instead
                "useMaterialParam"   : solverContext["useMaterialParameters"],
                "useMaterialParamVect" : solverContext["useMaterialParametersVect"],  #TODO JMG: legacy
                "useCERKGuess"       : solverContext["useCERK"],
                "useSplitCK"         : solverContext["useSplitCK"],
                "useGaussLobatto"    : solverContext["basis"] == "lobatto",
                "predictorRecompute" : solverContext["predictorRecompute"],
                "useVectPDE"         : solverContext["useVectPDE"],
                "useAoSoA2"          : solverContext["useAoSoA2"],
                # Optional int parameters (may set redundant flags)
                "usePointSources"    : solverContext["numberOfPointSources"] if solverContext["numberOfPointSources"] > 0 else -1,
                "tempVarsOnStack"    : solverContext["tempVarsOnStack"]
            }
        elif solverContext["kernelType"] == "limiter":
            kernelgeneratorContext = {
                # Mandatory parameters
                "kernelType"         : "limiter",
                "pathToApplication"  : solverContext["outputPath"],
                "pathToOptKernel"    : solverContext["optKernelPath"],
                "namespace"          : solverContext["optNamespace"],
                "solverName"         : solverContext["project"] + "::" + solverContext["solver"],
                "numberOfVariables"  : solverContext["numberOfVariables"],
                "numberOfParameters" : solverContext["numberOfMaterialParameters"],
                "order"              : solverContext["order"],
                "dimension"          : solverContext["dimensions"],
                "architecture"       : solverContext["architecture"],
                "numberOfObservable" : solverContext.get("numberOfDMPObservables", -1), #not set if not limiterSolver
                "ghostLayerWidth"    : solverContext.get("ghostLayerWidth", 0), #not set if not limiterSolver
                "limPatchSize"       : solverContext.get("patchSize",-1),
                # Optional bool parameters (may set redundant flags and default false flag)
                "useGaussLobatto"    : solverContext["basis"] == "lobatto",
                # Optional int parameters (may set redundant flags)
                "tempVarsOnStack"    : solverContext["tempVarsOnStack"]
            }
        elif solverContext["kernelType"] == "fv":
            kernelgeneratorContext = {
                # Mandatory parameters
                "kernelType"         : "fv",
                "pathToApplication"  : solverContext["outputPath"],
                "pathToOptKernel"    : solverContext["optKernelPath"],
                "namespace"          : solverContext["optNamespace"],
                "solverName"         : solverContext["project"] + "::" + solverContext["solver"],
                "numberOfVariables"  : solverContext["numberOfVariables"],
                "numberOfParameters" : solverContext["numberOfMaterialParameters"],
                "patchSize"          : solverContext["patchSize"],
                "dimension"          : solverContext["dimensions"],
                "finiteVolumesType"  : solverContext["finiteVolumesType"],
                "architecture"       : solverContext["architecture"],
                "slopeLimiter"       : solverContext["slopeLimiter"],
                # Optional bool parameters (may set redundant flags and default false flag)
                "useFlux"            : solverContext["useFlux"],
                "useViscousFlux"     : solverContext["useViscousFlux"],
                "useNCP"             : solverContext["useNCP"],
                "useSource"          : solverContext["useSource"],
                "useFusedSource"     : solverContext["useFusedSource"],
                "useMaterialParam"   : solverContext["useMaterialParameters"],
                "useRobustDiagLim"   : solverContext["useRobustDiagonalLimiting"],
                # Optional int parameters (may set redundant flags)
                "usePointSources"    : solverContext["numberOfPointSources"] if solverContext["numberOfPointSources"] > 0 else -1,
                "tempVarsOnStack"    : solverContext["tempVarsOnStack"]
            }
        else:
            raise ValueError("KernelType '"+context["kernelType"]+"' not supported by the kernelgenerator")
        # call the kernelgenerator with the given context
        kernelgeneratorController = kernelgenerator.Controller(kernelgeneratorContext)
        kernelgeneratorController.generateCode()
        
        # if verbose print the associated command line
        kernelgeneratorContext["commandLine"] = kernelgeneratorController.commandLine
        
        return solverContext["optKernelPath"], kernelgeneratorContext
