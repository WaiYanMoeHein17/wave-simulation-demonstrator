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
# Generate a ccph with getter to the parameters used by the kernel generator
#


from .abstractModelBaseClass import AbstractModelBaseClass


class ConfigurationParametersModel(AbstractModelBaseClass):

    def generateCode(self):
        
        # shortcut
        nVar     = self.context["nVar"]
        nVarPad  = self.context["nVarPad"]
        nPar     = self.context["nPar"]
        nParPad  = self.context["nParPad"]
        nData    = self.context["nData"]
        nDataPad = self.context["nDataPad"]
        nDof     = self.context["nDof"]
        nDofPad  = self.context["nDofPad"]
        nDof2Pad = self.context["nDof2Pad"]
        nDof3D    = self.context["nDof3D"]
        nDim     = self.context["nDim"]
        
        if self.context["kernelType"] == "aderdg":
            self.context["isLinearCText"]          = "true" if self.context["isLinear"]    else "false" #c++ true/false instead of True/False
            self.context["useFluxFactor"]          = 1 if self.context["useFlux"]          else 0
            self.context["useNCPFactor"]           = 1 if self.context["useNCP"]           else 0
            self.context["useSourceFactor"]        = 1 if self.context["useSource"]        else 0
            self.context["usePointSourcesFactor"]  = 1 if self.context["usePointSources"]  else 0
            self.context["useMaterialParamFactor"] = 1 if self.context["useMaterialParam"] else 0
            self.context["useSourceOrNCPFactor"]   = 1 if self.context["useSourceOrNCP"]   else 0
            
            # SPT buffer sizes
            # default value
            self.context["lQiSize"]   =   -1
            self.context["lQiNextSize"] = -1
            self.context["lPiSize"]   =   -1
            self.context["lFiSize"]   =   -1
            self.context["lSiSize"]   =   -1
            self.context["lQhiSize"]  =   -1
            self.context["lFhiSize"]  =   -1
            self.context["lShiSize"]  =   -1
            self.context["gradQSize"] =   -1
            self.context["PSiSize"]   =   -1
            self.context["PSiDerivativeSize"] =   -1
            if self.context["isLinear"]:
                if(self.context["useSplitCK"]):
                    if self.context["useVectPDE"]:
                        if self.context["useAoSoA2"]:
                            # Linear + split CK vect AoSoA2
                            self.context["lQiSize"]   = nDof3D*nVar*nDof2Pad
                            self.context["lQiNextSize"] = nDof3D*nVar*nDof2Pad
                            self.context["lPiSize"]   = nDof3D*nPar*nDof2Pad
                            self.context["lQhiSize"]  = nDof3D*nVar*nDof2Pad
                            self.context["lFhiSize"]  = nDof3D*nVar*nDof2Pad
                            self.context["gradQSize"] = nDof3D*nVar*nDof2Pad
                            self.context["PSiSize"]   = nDof*nDof3D*nVar*nDof2Pad
                        else:
                            # Linear + split CK vect (AoSoA)
                            self.context["lQiSize"]   = nDof3D*nDof*nVar*nDofPad
                            self.context["lQiNextSize"] = nDof3D*nDof*nVar*nDofPad
                            self.context["lPiSize"]   = nDof3D*nDof*nPar*nDofPad
                            self.context["lQhiSize"]  = nDof3D*nDof*nVar*nDofPad
                            self.context["lFhiSize"]  = nDof3D*nDof*nVar*nDofPad
                            self.context["gradQSize"] = nDof3D*nDof*nVar*nDofPad
                            self.context["PSiSize"]   = nDof*nDof3D*nDof*nVar*nDofPad
                    else:
                        # Linear + split CK scalar
                        self.context["lQiSize"]   = nVarPad*(nDof**nDim)
                        self.context["lQiNextSize"] = nVarPad*(nDof**nDim)
                        self.context["lPiSize"]   = nParPad*(nDof**nDim)
                        self.context["lQhiSize"]  = nVarPad*(nDof**nDim)
                        self.context["lFhiSize"]  = nVarPad*(nDof**nDim)
                        self.context["gradQSize"] = nVarPad*(nDof**nDim)
                        self.context["PSiSize"]   = nDof*(nDof**nDim)*nVarPad
                else: 
                    # default linear
                    self.context["lQiSize"]   = nDataPad*(nDof**nDim)*(1+nDof)
                    self.context["lQhiSize"]  = nDataPad*(nDof**nDim)
                    self.context["lFiSize"]   = nDim*nDof*(nDof**nDim)*nVarPad#nVarPad*(nDof**(nDim+1))*(2*nDim+1) # Todo JMG see if 2dim+1 or 2dim
                    self.context["lFhiSize"]  = nVarPad*(nDof**nDim)*nDim
                    if self.context["useSource"]:
                        self.context["lSiSize"]   = nVarPad*(nDof**(nDim+1))
                        self.context["lShiSize"]  = nVarPad*(nDof**nDim)
                    if self.context["useNCP"]:
                        self.context["gradQSize"] = nVarPad*(nDof**nDim)*nDim
                    if self.context["usePointSources"]:
                        self.context["PSiSize"]   = (nDof+1)*(nDof**nDim)*nVarPad
                        self.context["PSiDerivativeSize"] = self.context["PSiSize"]
            else:
            # nonlinear
                if self.context["predictorRecompute"]:
                    if self.context["useVectPDE"]:
                        self.context["lQiSize"]   = nDofPad*nVar*(nDof**nDim)
                        self.context["lQhiSize"]  = nDofPad*nVar*nDof*nDof3D
                        if nPar > 0:
                            self.context["lPiSize"]   = nDofPad*nPar*nDof*nDof3D
                        if self.context["useFlux"]:
                            self.context["lFhiSize"]  = nDofPad*nVar*nDof*nDof3D*nDim
                        if self.context["useSource"] or self.context["useNCP"]:
                            self.context["lShiSize"]  = nDofPad*nVar*nDof*nDof3D
                        if self.context["useNCP"] or self.context["useViscousFlux"]:
                            self.context["gradQSize"] = nDofPad*nVar*nDof*nDof3D*nDim
                    else: #scalar predictorRecompute
                        self.context["lQiSize"]   = nVarPad*(nDof**(nDim+1))
                        self.context["lQhiSize"]  = nVarPad*(nDof**nDim)
                        if nPar > 0:
                            self.context["lPiSize"]   = nParPad*(nDof**nDim)
                        if self.context["useFlux"]:
                            self.context["lFhiSize"]  = nVarPad*(nDof**nDim)*nDim
                        if self.context["useSource"] or self.context["useNCP"]:
                            self.context["lShiSize"]  = nVarPad*(nDof**nDim)
                        if self.context["useNCP"] or self.context["useViscousFlux"]:
                            self.context["gradQSize"] = nVarPad*(nDof**nDim)*nDim
                else: # default nonlinear
                    self.context["lQiSize"]   = nDataPad*(nDof**(nDim+1))
                    self.context["lQhiSize"]  = nDataPad*(nDof**nDim)
                    if self.context["useFlux"]:
                        self.context["lFiSize"]   = nVarPad*(nDof**(nDim+1))*nDim
                        self.context["lFhiSize"]  = nVarPad*(nDof**nDim)*nDim
                    if self.context["useSource"] or self.context["useNCP"]:
                        self.context["lSiSize"]   = nVarPad*(nDof**(nDim+1))
                        self.context["lShiSize"]  = nVarPad*(nDof**nDim)
                    if self.context["useNCP"] or self.context["useViscousFlux"]:
                        self.context["gradQSize"] = nVarPad*(nDof**nDim)*nDim
            
            # Face buffer size (Riemann)
            self.context["BndFaceSize"]      = nDataPad*(nDof*nDof3D)
            self.context["BndFaceTotalSize"] = 2*nDim*self.context["BndFaceSize"]
            self.context["BndFluxSize"]      = nVarPad*(nDof*nDof3D)
            self.context["BndFluxTotalSize"] = 2*nDim*self.context["BndFluxSize"]
        
        elif self.context["kernelType"] == "fv":
            # shortcut
            ghostLayerWidth = self.context["ghostLayerWidth"]
            
            # Ghostlayer bnd buffer size (BoundaryCondition)
            self.context["luhbndSize"] = nData*nDof*nDof3D*ghostLayerWidth
        
        self.render((self.context["kernelType"], "configurationParameters_cpph.template"), "ConfigurationParameters.cpph")
