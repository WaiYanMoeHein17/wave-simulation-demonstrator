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
#ifndef LIMITERPROJECTIONMATRICES_H_
#define LIMITERPROJECTIONMATRICES_H_

#include <set>

//TODO remove when generated
#include "kernels/KernelUtils.h"
#include <stdlib.h>

#include "tarch/Assertions.h"
#include "GaussLegendreBasis.h"
#include "GaussLobattoBasis.h"

namespace kernels {

// TODO Dominic: Old functions still needed?

///**
// * Initialises the lookup tables \p luh2lim, \p luh2lob
// * and \p lim2luh for the specified \p orders.
// *
// * \todo default implementation!
// *
// * \see freeLimiterProjectionMatrices
// */
//void initLimiterProjectionMatrices(const std::set<int>& orders);
//
///**
// * Frees the memory that was allocated for the lookup tables \p luh2lim, \p luh2lob
// * and \p lim2luh for the specified \p orders.
// *
// * \todo default implementation!
// *
// * \see initLimiterProjectionMatrices
// */
//void freeLimiterProjectionMatrices(const std::set<int>& orders);

extern double** uh2lim;
extern double** uh2lob;
extern double** lim2uh;

//TODO JMG remove when generated
void BaseFunc1D(double* phi, double xi, const int basisSize);
double* matrixInverse(int n, double* a);

// NEW

template <int n>
void matrixInverse(double (&ia) [n*n],const double (&a) [n*n]) {
  //TODO JMG remove when generated value
  std::fill_n(std::begin(ia),n*n,0.0);
  double c[n*n*2] = {0.0};

  idx2 idx(n,n);
  idx2 idxC(n,2*n);

  for(int i=0; i<n; i++) {
    for(int j=0; j<n; j++) {
      c[idxC(i,j)] = a[idx(j,i)];
    }
    for(int j; j<2*n; j++) {
      c[idxC(i,j)] = 0.;
    }
    c[idxC(i,i+n)] = 1.;
  }


  double tmp, piv, mlV;

  //Forward elimination and row swapping (if necessary)
  for(int i=0; i<n; i++) {
    int ml = i;
    mlV = std::abs(static_cast<double>(c[idxC(i,i)]));
    for(int j=i+1; j<n; j++) {
      if(std::abs(c[idxC(j,i)]) > mlV) {
        ml = j;
        mlV = c[idxC(j,i)];
      }
    }

    for(int k=0; k<2*n; k++) {
      tmp = c[idxC(ml,k)];
      c[idxC(ml,k)] = c[idxC(i,k)];
      c[idxC(i,k)] = tmp;
    }
    if(c[idxC(i,i)] == 0) {
      //logError("matrixInverse()", "Matrix is singular" );
    }
    piv = 1. / c[idxC(i,i)];
    for(int k=0; k<2*n; k++) {
      c[idxC(i,k)] *= piv;
    }
    for(int j=i+1; j<n; j++) {
      tmp = c[idxC(j,i)];
      for(int k=0; k<2*n; k++) {
        c[idxC(j,k)] -= tmp*c[idxC(i,k)];
      }
    }
  }

  //Back substitution
  for(int i=n-1; i>=0; i--) {
    for(int j=i-1; j>=0; j--) {
      tmp = c[idxC(j,i)];
      for(int k=0; k<2*n; k++) {
        c[idxC(j,k)] -= tmp*c[idxC(i,k)];
      }
    }
  }

  for(int i=0; i<n; i++) {
    for(int j=0; j<n; j++) {
      ia[idx(j,i)] = c[idxC(i,j+n)];
    }
  }
}

/**
 * For order \p order, evaluate all \p order+1 Gauss-Legendre basis functions.
 *
 * @param[inout] phi the values of the basis functions at point xi
 * @param[inout] xi the evaluation point.
 */
template <int order>
void gaussLegendreBasisFunctions(double (&phi) [order+1], double xi) {
  constexpr int basisSize = order+1;

  for(int i=0; i<basisSize; i++) {
    phi[i] = 1;
  }
  for(int m=0; m<basisSize; m++) {
    for(int j=0; j<basisSize; j++) {
      if(j == m)
        continue;
      else
        phi[m] = phi[m]*(xi - kernels::legendre::nodes[order][j])/(kernels::legendre::nodes[order][m]-kernels::legendre::nodes[order][j]);
    }
  }
}


/**
 * Compute the DG to FV projectors.
 * We simply evaluate the DG basis functions at the
 * centres of the finite volumes.
 *
 * @tparam order the polynomial order.
 * @tparam patchSize the size of the FV patch.
 *
 * @param[inout] dg2fv the DG to FV projector.
 */
template <int order,int patchSize>
void computeDG2FVProjector(double (&dg2fv) [(order+1)*patchSize]) {
  constexpr int basisSize    = order+1;
  constexpr int basisSizeLim = patchSize;
  constexpr double dxi = 1. / basisSizeLim;

  assertion2(basisSize < basisSizeLim,basisSize,basisSizeLim);

  double phi[basisSize] = {0.0};
  {
    idx2 idx(basisSize, basisSizeLim);
    for(int i=0; i<basisSizeLim; i++) {
      const double xLeft = i*dxi;
      for(int j=0; j<basisSize; j++) {
        const double xi = xLeft + dxi*kernels::legendre::nodes[order][j];
        gaussLegendreBasisFunctions<order>(phi, xi);
        for(int k=0; k<basisSize; k++) { //
          dg2fv[idx(k,i)] += kernels::legendre::weights[order][j] * phi[k];
        }
      }
    }
  }
}

/**
 * There are typically more finite volumes than DG basis functions per patch and per space dimension.
 * Therefore, we need to solve a minimisation problem to derive a DG solution.
 * We use a simple least-squares minisation subject to mass conservation, i.e.
 * the average DG solution over the patch must equal the average FV solution over the patch.
 *
 * @precondition The DG-to-FV projector must be computed before calling this function.
 *
 * @tparam order the polynomial order.
 * @tparam patchSize the size of the FV patch.
 *
 * @param[in]    dg2fv the DG to FV projector.
 * @param[inout] fv2dg the resulting FV to DG projector.
 *
 * References:
 * -----------
 *
 * Main reference:
 *
 * M. Dumbser et al., A posteriori subcell limiting of the discontinuous Galerkin finite element method for hyperbolic conservation laws,
 * Journal of Computational Physics, vol. 278, pp. 47–75, Dec. 2014.
 *
 * Constrained least-squares reference:
 *
 * M. Dumbser and M. Käser, Arbitrary high order non-oscillatory finite volume schemes on unstructured meshes for linear hyperbolic systems,
 * Journal of Computational Physics, vol. 221, no. 2, pp. 693–723, Feb. 2007.
 */
template <int order,int patchSize>
void computeFV2DGProjector(double (&fv2dg) [(order+1)*patchSize],const double (&dg2fv) [(order+1)*patchSize]) {
  constexpr int basisSize    = order+1;
  constexpr int basisSizeLim = patchSize;

  assertion2(basisSize < basisSizeLim,basisSize,basisSizeLim);

  double lsqm  [(basisSize+1)*(basisSize+1)] = {0.0};
  double lsqrhs[basisSizeLim*(basisSize+1) ] = {0.0};

  idx2 idxFV2DG(basisSizeLim, basisSize);
  idx2 idxLSQM((basisSize+1),(basisSize+1));
  idx2 idxLSQrhs(basisSizeLim,(basisSize+1));
  idx2 idxDG2FV(basisSize, basisSizeLim);
  const double dxi = 1.0 / basisSizeLim;
  for(int i=0; i<basisSize; i++) {
    for(int j=0; j<basisSize; j++) {
      lsqm[idxLSQM(i,j)] = 0.;
      for(int k=0; k<basisSizeLim; k++) {
        lsqm[idxLSQM(i,j)] += 2* dg2fv[idxDG2FV(i,k)] * dg2fv[idxDG2FV(j,k)];
      }
    }
    lsqm[idxLSQM(i,basisSize)] = kernels::legendre::weights[basisSize-1][i];
  }
  for(int i=0; i<basisSize; i++) {
    lsqm[idxLSQM(basisSize,i)] = -kernels::legendre::weights[basisSize-1][i];
  }
  lsqm[idxLSQM(basisSize,basisSize)] = 0.;

  double* ilsqm = matrixInverse(basisSize+1, lsqm);

  for(int i=0; i<basisSizeLim; i++) {
    for(int j=0; j<basisSize; j++) {
      lsqrhs[idxLSQrhs(i,j)] = 2*dg2fv[idxDG2FV(j,i)];
    }
    lsqrhs[idxLSQrhs(i,basisSize)] = dxi;
  }

  for(int i=0; i<basisSizeLim; i++) {
    for(int j=0; j<basisSize; j++) {
      fv2dg[idxFV2DG(i,j)] = 0.;
      for(int k=0; k<basisSize+1; k++) {
        fv2dg[idxFV2DG(i,j)] += ilsqm[idxLSQM(k,j)] * lsqrhs[idxLSQrhs(i,k)];
      }
    }
  }
}

/**
 * Compute a projector to switch from Gauss-Legendre basis functions
 * to Gauss-Lobatto basis functions.
 *
 * @tparam order the polynomial order.
 */
template <int order>
void computeLegendre2LobattoProjector(double (&leg2log) [(order+1)*(order+1)]) {
  constexpr int basisSize = order+1;

  double phi[basisSize] = {0.0};
  idx2 idx(basisSize, basisSize);

  for(int i=0; i<basisSize; i++) {
    gaussLegendreBasisFunctions<order>(phi, kernels::lobatto::nodes[order][i]);
    for(int j=0; j<basisSize; j++) {
      leg2log[idx(j,i)] = phi[j]; //Fortran: uh2lob(ii,:) = phi(:)
    }
  }
}

}

#endif //LIMITERPROJECTIONMATRICES_H_
