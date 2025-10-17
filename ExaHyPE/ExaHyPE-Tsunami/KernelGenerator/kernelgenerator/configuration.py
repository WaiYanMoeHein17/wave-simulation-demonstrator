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
# Static configuration parameters for the module
#
# @note
# requires python3

import sys
import os

class Configuration:

    ######################################
    ###### Configuration parameters ######
    ######################################
    # Change them if required

    # path to the root of ExaHyPe from this file
    pathToExaHyPERoot          = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

    # path to the gemm generator from this file
    pathToLibxsmmGemmGenerator = os.path.abspath(os.path.join(pathToExaHyPERoot, "Submodules", "libxsmm", "bin", "libxsmm_gemm_generator"))

    # path to eigen, will be symlinked into the kernels directory if eigen is used (see matmulLib)
    pathToEigen                = os.path.abspath(os.path.join(pathToExaHyPERoot, "Submodules", "eigen"))

    # path to jinja2
    pathToJinja2               = os.path.abspath(os.path.join(pathToExaHyPERoot, "Submodules", "jinja", "src"))

    # path to markupsafe
    pathToMarkupsafe           = os.path.abspath(os.path.join(pathToExaHyPERoot, "Submodules", "markupsafe", "src"))
    
    # internal path to templates
    pathToTemplate             = os.path.abspath(os.path.join(os.path.dirname(__file__), "templates"))

    # simd size of the accepted architectures
    simdWidth = { "noarch" : 1,
                  "wsm"    : 2,
                  "snb"    : 4,
                  "hsw"    : 4,
                  "knc"    : 8,
                  "knl"    : 8,
                  "skx"    : 8
                }

    # choose the BLAS library for the matmul: "None" (= C++ loops), "Libxsmm" or "Eigen"
    matmulLib = "Libxsmm"
    #matmulLib = "Eigen"
    #matmulLib = "None"

    # set to true to print models runtime
    runtimeDebug = False
    
    # prefetching settings
    # Experimental, not supported by all kernel
    # Will use prefetching to optimize tensor operation (prefetch the next slice of an LoG)
    prefetching = "None" # "None", "Inputs", "Outputs", "All"
    prefetchLevel = "_MM_HINT_T0" # intrisic _mm_prefetch locality hint (_MM_HINT_T0 = all level of cache), see compiler header xmmintrin.h
    cachelineSize = {
                  "noarch" : 8,
                  "wsm"    : 8,
                  "snb"    : 8,
                  "hsw"    : 8,
                  "knc"    : 8,
                  "knl"    : 8,
                  "skx"    : 8
    } # size for double, CPUs usuallly have 64B L1 Cache line



def checkPythonVersion():
    """check version. Python 3.6 required"""
    requiredVersion = (3,6)
    currentVersion  = sys.version_info
    if(requiredVersion > currentVersion):
        sys.exit("Requires Python 3.6 or newer. Abort.")



def checkDependencies():
    """check all dependencies are reachable from the configuration path"""
    # Check jinja
    sys.path.insert(1, Configuration.pathToJinja2)
    sys.path.insert(1, Configuration.pathToMarkupsafe)
    import jinja2
    # Remove added path
    sys.path.remove(Configuration.pathToJinja2)
    sys.path.remove(Configuration.pathToMarkupsafe)
