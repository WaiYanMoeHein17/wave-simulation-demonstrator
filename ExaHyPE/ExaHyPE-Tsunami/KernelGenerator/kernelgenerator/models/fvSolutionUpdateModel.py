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
# Generates the solutionUpdate Kernel
#


from .abstractModelBaseClass import AbstractModelBaseClass


class FVSolutionUpdateModel(AbstractModelBaseClass):

    def generateCode(self):
        self.context["nDofGS"] = self.context["nDof"]+2*(self.context["ghostLayerWidth"]-1); # patchsize + (ghostlayers - outermost ghostlayer)
        self.context["nDofGS3D"] = 1 if self.context["nDim"] == 2 else self.context["nDofGS"]
        self.render(("fv", "solutionUpdate_cpp.template"), "solutionUpdate.cpp")
