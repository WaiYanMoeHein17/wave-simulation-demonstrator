import os
import sys
import copy

from .models import *
from .toolkitHelper import ToolkitHelper

class SolverController:
    def __init__(self, solverSpec, baseContext):
        self.solverSpec = solverSpec
        self.baseContext = baseContext

        # ADER-DG CFL factors        
        # Up to order 4, computed via von Neumann analysis [1]; above, determined empirically by M. Dumbser)
        # [1] M. Dumbser, D. S. Balsara, E. F. Toro, and C.-D. Munz, 
        # ‘A unified framework for the construction of one-step finite volume and discontinuous Galerkin schemes on unstructured meshes’, Journal of Computational Physics, vol. 227, no. 18, pp. 8209–8253, Sep. 2008.
        self.cflADER = [1.0,   0.33,  0.17, 0.1,  0.069, 0.045, 0.038, 0.03, 0.02, 0.015]
        # Values for orders > 9 are found experimentally, too.
        p = 10
        #for cflCorrection in [0.1995, 0.13965,  0.097755, 0.0684285, 0.04789995, 0.033529965]:
        for cflCorrection in [0.18525, 0.1204125, 0.078268125, 0.05087428125, 0.0330682828125, 0.021494383828125]:
            self.cflADER.append(cflCorrection / (2*p+1))
            p += 1
        self.cflCorrectionADER = self.cflADER.copy()
        for p,cfl in enumerate(self.cflADER):
            self.cflCorrectionADER[p] = cfl * (2*p+1);
        self.cflSafetyFactor = 0.9


    def processModelOutput(self, output, contextsList, logger):
        """Standard model output is (paths, context)
        
        Log the path of the generated file (== not None) using the logger and
        store the context used in the contextList
        
        return the context
        """
        paths, context = output
        for path in filter(None, paths):
            logger.info("Generated '"+path+"'")
        if "kernelgeneratorContext" in context:
            logger.info("Kernelgenerator used, command line to get the same result: "+context["kernelgeneratorContext"]["commandLine"])
        contextsList.append(context)
        
        return context


    def run(self, logger):
        solverContextsList = []
        for i,solver in enumerate(self.solverSpec):
            logger.info("Generating solver[%d] = %s..." % (i, solver["name"]))
            if solver["type"]=="ADER-DG": 
                model = solverModel.SolverModel(self.buildADERDGSolverContext(solver,logger))
                solverContext = self.processModelOutput(model.generateCode(), solverContextsList, logger)
            elif solver["type"]=="Finite-Volumes":
                context = self.buildFVSolverContext(solver)
                if type(context["patchSize"]) != int: # TODO move to proper validate
                    logger.error("{}: 'patch_size' must be an integer; it is '{}'.".format(context["solver"],context["patchSize"]))
                    raise ValueError("'patch_size' must be an integer")
                model = solverModel.SolverModel(context)
                solverContext = self.processModelOutput(model.generateCode(), solverContextsList, logger)
            elif solver["type"]=="Limiting-ADER-DG":
                aderdgContext = self.buildADERDGSolverContext(solver, logger)
                fvContext     = self.buildFVSolverContext(solver)
                context       = self.buildLimitingADERDGSolverContext(solver)
                # modifications
                fvContext["solver"]         = context["FVSolver"]
                fvContext["solverType"]     = "Finite-Volumes"
                fvContext["abstractSolver"] = context["FVAbstractSolver"]
                
                # patch size vs ADER-DG CFL
                order = aderdgContext["order"];
                if fvContext["patchSize"] == "default":
                    fvContext["patchSize"] = 2 * order + 1
                elif fvContext["patchSize"] == "max":
                    fvContext["patchSize"] = int( 1.0/self.cflCorrectionADER[order] * (2 * order + 1) )
                else:
                    # scale the ADER-DG method's CFL safety factor according to the chosen patch size
                    aderdgContext["CFL"] = aderdgContext["CFL"] * min(self.cflADER[order], 1.0/fvContext["patchSize"]) / self.cflADER[order];

                # forward informations
                context["patchSize"] = fvContext["patchSize"]
                context["ghostLayerWidth"] = fvContext["ghostLayerWidth"]

                aderdgContext["solver"]                 = context["ADERDGSolver"]
                aderdgContext["solverType"]             = "ADER-DG"
                aderdgContext["abstractSolver"]         = context["ADERDGAbstractSolver"]
                aderdgContext["numberOfDMPObservables"] = context["numberOfDMPObservables"]
                aderdgContext["ghostLayerWidth"]        = fvContext["ghostLayerWidth"]
                
                self.addKernelgeneratorPathAndNamespace(aderdgContext) # refresh path and namespace
                self.addKernelgeneratorPathAndNamespace(fvContext) # refresh path and namespace
                
                # generate all solver files
                model = solverModel.SolverModel(fvContext)
                context["fvContext"] = self.processModelOutput(model.generateCode(), [], logger) #don't register context
                model = solverModel.SolverModel(aderdgContext)
                context["aderdgContext"] = self.processModelOutput(model.generateCode(), [], logger) #don't register context
                if "kernelgeneratorContext" in context["aderdgContext"]:
                    context["kernelgeneratorContext"] = context["aderdgContext"]["kernelgeneratorContext"] #move kernelgencontext one up if it exists
                    # Add missing, TODO JMG make cleaner
                    context["basis"] = context["aderdgContext"]["basis"] 
                    context["tempVarsOnStack"] = context["aderdgContext"]["tempVarsOnStack"]
                model = solverModel.SolverModel(context)
                solverContext = self.processModelOutput(model.generateCode(), solverContextsList, logger)
                
            solverContext["plotters"] = []
            for j, plotter in enumerate(solver.get("plotters",[])):
                logger.info("Generating plotter[%d]= %s ,for solver[%d]= %s..." % (j, plotter["name"], i, solver["name"]))
                plotModel = plotterModel.PlotterModel(self.buildPlotterContext(solver,plotter))
                self.processModelOutput(plotModel.generateCode(), solverContext["plotters"], logger)
            
        logger.info("Solver generation... done")
        
        return solverContextsList

    def buildBaseSolverContext(self, solver):
        context = copy.copy(self.baseContext)
        
        context["solverType"]     = solver["type"]
        context["solver"]         = solver["name"]
        context["abstractSolver"] = "Abstract"+context["solver"]
        
        nVar          = ToolkitHelper.count_variables(ToolkitHelper.parse_variables(solver,"variables"))
        nParam        = ToolkitHelper.count_variables(ToolkitHelper.parse_variables(solver,"material_parameters"))
        nGlobalObs    = ToolkitHelper.count_variables(ToolkitHelper.parse_variables(solver,"global_observables"))
        nPointSources = solver["point_sources"] if type(solver.get("point_sources",[])) is int else len(solver.get("point_sources",[]))
        
        context["numberOfVariables"]          = nVar
        context["numberOfParameters"]         = nParam
        context["numberOfMaterialParameters"] = nParam
        context["numberOfGlobalObservables"]  = nGlobalObs
        context["numberOfPointSources"]       = nPointSources
        context["CFL"]                        = solver["cfl"] # default value is 0.9

        # variables access class
        context["variablesMap"]  = ToolkitHelper.parse_variables(solver,"variables")
        if nParam>0:
            parametersMap = ToolkitHelper.parse_variables(solver,"material_parameters")
            # Increase offset of parameters, they are located directly after variables.
            def increaseOffset(param):
                param["offset"] += nVar
                return param
            parametersMap = [increaseOffset(param) for param in parametersMap]
            context["variablesMap"] += parametersMap
        context["variablesMapSize"] = len(context["variablesMap"])
        context["variables_as_str"] = ToolkitHelper.variables_to_str(solver,"variables")
        context["material_parameters_as_str"]  = ToolkitHelper.variables_to_str(solver,"material_parameters")
        
        context["globalObservablesMap"]      = ToolkitHelper.parse_variables(solver,"global_observables")
        context["global_observables_as_str"] = ToolkitHelper.variables_to_str(solver,"global_observables")
        
        context["range_0_nVar"]          = range(0,nVar)
        context["range_0_nVarParam"]     = range(0,nVar+nParam)
        context["range_0_nGlobalObs"]    = range(0,nGlobalObs)    # nGlobalObs might be 0
        context["range_0_nPointSources"] = range(0,nPointSources) # nPointSources might be 0
        
        context["namingSchemes"]={} # TODO read from spec
        
        return context


    def buildADERDGSolverContext(self, solver, logger):
        context = self.buildBaseSolverContext(solver)
        context["type"] = "ADER-DG"
        context.update(self.buildADERDGKernelContext(solver["aderdg_kernel"]))
        context.update(self.buildKernelOptimizationContext(solver["aderdg_kernel"], context))
        self.addKernelgeneratorPathAndNamespace(context)
        
        context["order"]                  = solver["order"]
        if int(context["order"]) > 9:
            logger.warning("Support for orders greater than 9 is currently only experimental!")

        context["PNPM"]                   = self.cflADER[int(solver["order"])]
        context["numberOfDMPObservables"] = -1 # overwrite if called from LimitingADERDGSolver creation

        return context


    def buildLimitingADERDGSolverContext(self, solver):
        context = self.buildBaseSolverContext(solver)
        context["type"]                   = "Limiting-ADER-DG"
        context["order"]                  = solver["order"]
        context["numberOfDMPObservables"] = solver["limiter"]["dmp_observables"]
        context["implementation"]         = solver["limiter"].get("implementation","generic")
        context["ADERDGSolver"]           = solver["name"]+"_ADERDG"
        context["FVSolver"]               = solver["name"]+"_FV"
        context["ADERDGAbstractSolver"]   = "Abstract"+solver["name"]+"_ADERDG"
        context["FVAbstractSolver"]       = "Abstract"+solver["name"]+"_FV"
        self.addKernelgeneratorPathAndNamespace(context)
        
        return context


    def buildFVSolverContext(self, solver):
        context = self.buildBaseSolverContext(solver)
        context["type"] = "Finite-Volumes"
        context.update(self.buildFVKernelContext(solver["fv_kernel"]))
        context.update(self.buildKernelOptimizationContext(solver["fv_kernel"], context))
        context["patchSize"] = solver.get("patch_size","default") # overwrite if called from LimitingADERDGSolver creation
        self.addKernelgeneratorPathAndNamespace(context)
        
        return context


    def buildADERDGKernelContext(self, kernel):
        context = {}
        context["implementation"]          = kernel.get("implementation","generic")
        context["useMaxPicardIterations"]  = kernel.get("space_time_predictor",{}).get("fix_picard_iterations",False)!=False
        context["tempVarsOnStack"]         = kernel.get("allocate_temporary_arrays","heap")=="stack" 
        context["patchwiseAdjust"]         = kernel.get("adjust_solution","pointwise")=="patchwise" 
        context["transformRiemannData"]    = kernel.get("transform_riemann_data",False)==True
        context["language"]                = kernel.get("language","C").lower()
        context["basis"]                   = kernel.get("basis","Legendre").lower()
        context["isLinear"]                = not kernel.get("nonlinear",True)
        context["isLinear_s"]              = "true" if not kernel.get("nonlinear",True) else "false"
        context["isNonlinear"]             = kernel.get("nonlinear",True)
        context["linearOrNonlinear"]       = "Linear" if context["isLinear"] else "Nonlinear"
        context["isFortran"]               = kernel.get("language",False)=="Fortran" 
        context["useCERK"]                 = kernel.get("space_time_predictor",{}).get("cerkguess",False)
        context["useSplitCK"]              = kernel.get("space_time_predictor",{}).get("split_ck",False)
        context["noTimeAveraging"]         = kernel.get("space_time_predictor",{}).get("notimeavg",False)
        context["noTimeAveraging_s"]       = "true" if kernel.get("space_time_predictor",{}).get("notimeavg",False) else "false"
        context["predictorRecompute"]      = kernel.get("space_time_predictor",{}).get("predictor_recompute",False)
        context["useVectPDE"]              = kernel.get("space_time_predictor",{}).get("vectorise_terms",False)
        context["useAoSoA2"]               = kernel.get("space_time_predictor",{}).get("AoSoA2_layout",False)
        context.update(self.buildKernelTermsContext(kernel["terms"]))
        return context


    def buildFVKernelContext(self,kernel):
        context = {}
        ghostLayerWidth = { "godunov" : 1, "musclhancock" : 2 }
        context["finiteVolumesType"]           = kernel["scheme"].replace("robust","")
        context["ghostLayerWidth"]             = ghostLayerWidth[context["finiteVolumesType"]]
        context["useRobustDiagonalLimiting"] = ("robust" in kernel["scheme"])
        context["useRobustDiagonalLimiting_s"] = "true" if context["useRobustDiagonalLimiting"] else "false"
        context["slopeLimiter"]   = kernel.get("slope_limiter","minmod")

        context["implementation"]  = kernel.get("implementation","generic")
        context["tempVarsOnStack"] = kernel.get("allocate_temporary_arrays","heap")=="stack" 
        context["patchwiseAdjust"] = kernel.get("adjust_solution","pointwise")=="patchwise" 
        context.update(self.buildKernelTermsContext(kernel["terms"]))
        return context


    def buildKernelTermsContext(self,terms):
        context = {}
        context["useFlux"]                  = "flux" in terms or "viscous_flux" in terms
        context["useFlux_s"]                = "true" if context["useFlux"] else "false"
        context["useSource"]                = "source" in terms
        context["useSource_s"]              = "true" if context["useSource"] else "false"
        context["useNCP"]                   = "ncp" in terms
        context["useNCP_s"]                 = "true" if context["useNCP"] else "false"
        context["usePointSources"]          = "point_sources" in terms
        context["usePointSources_s"]        = "true" if context["usePointSources"] else "false"
        context["useMaterialParameters"]    = "material_parameters" in terms
        context["useMaterialParameters_s"]  = "true" if context["useMaterialParameters"] else "false"
        context["useViscousFlux"]           = "viscous_flux" in terms
        context["useViscousFlux_s"]         = "true" if context["useViscousFlux"] else "false"
        
        return context

    def buildPlotterContext(self,solver,plotter):
        context = self.buildBaseSolverContext(solver)

        context["plotter"]          = plotter["name"]
        context["writtenUnknowns"]  = ToolkitHelper.count_variables(ToolkitHelper.parse_variables(plotter,"variables"))
        context["variables_as_str"] = ToolkitHelper.variables_to_str(plotter,"variables")
        context["plotterType"]      = plotter["type"] if type(plotter["type"]) is str else "::".join(plotter["type"])
        
        context["headerPath"]       = os.path.join(context["plotterSubDirectory"],(context["plotter"]+".h"))
        if context["plotterSubDirectory"] != "":
            context["outputPath"]   = os.path.join(context["outputPath"], context["plotterSubDirectory"])

            # This might be the wrong place, but we need to make sure the plotterSubDirectory exists.
            os.makedirs(context["outputPath"], exist_ok=True)
            
        return context


    def buildKernelOptimizationContext(self, kernel, solverContext):
        optimizations = kernel.get("optimised_kernel_debugging",[]) + kernel.get("optimised_terms",[])
        context = {}
        context["useConverter"]       = "converter"        in optimizations
        context["countFlops"]         = "flops"            in optimizations
        context["useFluxVect"]        = "flux_vect"        in optimizations
        context["useFusedSource"]     = "fusedsource"      in optimizations
        context["useFusedSourceVect"] = "fusedsource_vect" in optimizations
        context["useSourceVect"]      = "source_vect"      in optimizations
        context["useNCPVect"]         = "ncp_vect"         in optimizations
        context["useMaterialParametersVect"] = "material_parameters_vect" in optimizations
        
        return context
        
    def addKernelgeneratorPathAndNamespace(self, context):
        kernelType = "aderdg" if context["type"] == "ADER-DG" else ("limiter" if context["type"] == "Limiting-ADER-DG" else "fv")
        context["kernelType"]         = kernelType
        context["optKernelPath"]      = os.path.join("kernels", context["project"] + "_" + context["solver"], kernelType)
        context["optNamespace"]       = context["project"] + "::" + context["solver"] + "_kernels::"+kernelType
