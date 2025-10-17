/**
 * This file is part of the ExaHyPE project.
 * Copyright (c) 2016  http://exahype.eu
 * All rights reserved.
 *
 * The project has received funding from the European Union's Horizon 
 * 2020 research and innovation programme under grant agreement
 * No 671698. For copyrights and licensing, please consult the webpage.
 *
 * Released under the BSD 3 Open Source License.
 * For the full license text, see LICENSE.txt
 **/
#ifndef GAUSSLEGENDRE_H_
#define GAUSSLEGENDRE_H_

#include "kernels/KernelUtils.h"

// Snippet generated with ExaHyPE-Engine/Miscellaneous/aderdg/generateLookupTable.py

#include "kernels/generated/GaussLegendreBasis.hsnippet"

namespace kernels {
namespace legendre {

/**
 * Returns the exact (point-wise) interpoland for position for x
 *
 * If you interpolate onto a equidistant grid, you can alternatively use the
 * array equidistantGridProjector1d that holds all mapping quantities in a
 * precomputed table and thus is faster.
 *
 * @param offsetOfPatch Array of doubles of size DIMENSIONS. If you use
 *          Peano's tarch::la::Vector class, apply its data() function to
 *          get hold of a pointer to be passed into this function
 * @param x Position where to evaluate the function. This has to be a point
 *          within the patch
 * @param numberOfUnknowns Number of unknowns held per grid point withint the
 *          patch
 * @param unknown Which unknown to evaluate
 * @param order   Which order is used
 * @param u       Pointer to unknowns
 */
double interpolate(
    const double* offsetOfPatch,
    const double* sizeOfPatch,
    const double* x,
    int           numberOfUnknowns,
    int           unknown,
    int           order,
    const double* u
);

}
}

#endif
