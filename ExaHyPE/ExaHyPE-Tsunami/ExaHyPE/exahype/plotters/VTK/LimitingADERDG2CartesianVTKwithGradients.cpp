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
 
#include "LimitingADERDG2CartesianVTKwithGradients.h"
#include "tarch/parallel/Node.h"

// @todo 16/05/03:Dominic Etienne Charreir Plotter depends now on kernels.
// Should thus be placed in kernel module or the solver
// should provide a function that computes solution values
// at equidistant grid points
#include "kernels/GaussLegendreBasis.h"
#include "peano/utils/Loop.h"


#include "tarch/plotter/griddata/unstructured/vtk/VTKTextFileWriter.h"
#include "tarch/plotter/griddata/unstructured/vtk/VTKBinaryFileWriter.h"
#include "tarch/plotter/griddata/unstructured/vtk/VTUTextFileWriter.h"
#include "tarch/plotter/griddata/unstructured/vtk/VTUBinaryFileWriter.h"

#include "exahype/solvers/LimitingADERDGSolver.h"


// VTK subclasses
std::string exahype::plotters::LimitingADERDG2CartesianVerticesVTKAsciiwithGradients::getIdentifier() {
  return "vtk::Cartesian::vertices::limited::ascii::wG";
}
exahype::plotters::LimitingADERDG2CartesianVerticesVTKAsciiwithGradients::LimitingADERDG2CartesianVerticesVTKAsciiwithGradients(
    exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing,
    const int ghostLayerWidth):
  LimitingADERDG2CartesianVTKwithGradients(postProcessing,ghostLayerWidth,PlotterType::ASCIIVTK,false) {
}

std::string exahype::plotters::LimitingADERDG2CartesianVerticesVTKBinarywithGradients::getIdentifier() {
  return "vtk::Cartesian::vertices::limited::binary::wG";
}
exahype::plotters::LimitingADERDG2CartesianVerticesVTKBinarywithGradients::LimitingADERDG2CartesianVerticesVTKBinarywithGradients(
    exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing,
    const int ghostLayerWidth):
    LimitingADERDG2CartesianVTKwithGradients(postProcessing,ghostLayerWidth,PlotterType::BinaryVTK,false) {
}

std::string exahype::plotters::LimitingADERDG2CartesianCellsVTKAsciiwithGradients::getIdentifier() {
  return "vtk::Cartesian::cells::limited::ascii::wG";
}
exahype::plotters::LimitingADERDG2CartesianCellsVTKAsciiwithGradients::LimitingADERDG2CartesianCellsVTKAsciiwithGradients(
    exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing,
    const int ghostLayerWidth):
    LimitingADERDG2CartesianVTKwithGradients(postProcessing,ghostLayerWidth,PlotterType::ASCIIVTK,true) {
}

std::string exahype::plotters::LimitingADERDG2CartesianCellsVTKBinarywithGradients::getIdentifier() {
 return "vtk::Cartesian::cells::limited::binary::wG";
}
exahype::plotters::LimitingADERDG2CartesianCellsVTKBinarywithGradients::LimitingADERDG2CartesianCellsVTKBinarywithGradients(
    exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing,
    const int ghostLayerWidth):
    LimitingADERDG2CartesianVTKwithGradients(postProcessing,ghostLayerWidth,PlotterType::BinaryVTK,true) {
}

// VTU subclasses
std::string exahype::plotters::LimitingADERDG2CartesianVerticesVTUAsciiwithGradients::getIdentifier() {
  return "vtu::Cartesian::vertices::limited::ascii::wG";
}
exahype::plotters::LimitingADERDG2CartesianVerticesVTUAsciiwithGradients::LimitingADERDG2CartesianVerticesVTUAsciiwithGradients(
    exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing,
    const int ghostLayerWidth):
    LimitingADERDG2CartesianVTKwithGradients(postProcessing,ghostLayerWidth,PlotterType::ASCIIVTU,false) {
}

std::string exahype::plotters::LimitingADERDG2CartesianVerticesVTUBinarywithGradients::getIdentifier() {
  return "vtu::Cartesian::vertices::limited::binary::wG";
}
exahype::plotters::LimitingADERDG2CartesianVerticesVTUBinarywithGradients::LimitingADERDG2CartesianVerticesVTUBinarywithGradients(
    exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing,
    const int ghostLayerWidth):
    LimitingADERDG2CartesianVTKwithGradients(postProcessing,ghostLayerWidth,PlotterType::BinaryVTU,false) {
}

std::string exahype::plotters::LimitingADERDG2CartesianCellsVTUAsciiwithGradients::getIdentifier() {
  return "vtu::Cartesian::cells::limited::ascii::wG";
}
exahype::plotters::LimitingADERDG2CartesianCellsVTUAsciiwithGradients::LimitingADERDG2CartesianCellsVTUAsciiwithGradients(
    exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing,
    const int ghostLayerWidth):
    LimitingADERDG2CartesianVTKwithGradients(postProcessing,ghostLayerWidth,PlotterType::ASCIIVTU,true) {
}

std::string exahype::plotters::LimitingADERDG2CartesianCellsVTUBinarywithGradients::getIdentifier() {
 return "vtu::Cartesian::cells::limited::binary::wG";
}
exahype::plotters::LimitingADERDG2CartesianCellsVTUBinarywithGradients::LimitingADERDG2CartesianCellsVTUBinarywithGradients(
    exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing,
    const int ghostLayerWidth):
    LimitingADERDG2CartesianVTKwithGradients(postProcessing,ghostLayerWidth,PlotterType::BinaryVTU,true) {
}

tarch::logging::Log exahype::plotters::LimitingADERDG2CartesianVTKwithGradients::_log("exahype::plotters::LimitingADERDG2CartesianVTKwithGradients");

exahype::plotters::LimitingADERDG2CartesianVTKwithGradients::LimitingADERDG2CartesianVTKwithGradients(
    exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing,
    const int ghostLayerWidth,
    PlotterType plotterType,
    const bool plotCells)
  :
  Device(postProcessing),
  _plotterType(plotterType),
  _plotCells(plotCells),
  _ghostLayerWidth(ghostLayerWidth)
{
  if ( !plotCells ) {
    logError("LimitingADERDG2CartesianVTKwithGradients(...)", "Currently, there exists only an implementation of this plotter" <<
        "that writes cell data. Writing vertex data is currently not supported. " <<
        "You can increase the resolution of the plotted ADER-DG subcell averages via the plotter parameter " <<
        "'resolution'.");
    std::terminate();
  }
}


void exahype::plotters::LimitingADERDG2CartesianVTKwithGradients::init(
  const std::string& filename,
  int                orderPlusOne,
  int                unknowns,
  int                writtenUnknowns,
  exahype::parser::ParserView plotterParameters
) {
  _filename          = filename;
  _order             = orderPlusOne-1;
  _solverUnknowns    = unknowns;
  _plotterParameters = plotterParameters;
  _patchWriter       = nullptr;
  _writtenUnknowns   = writtenUnknowns;

  _slicer = Slicer::bestFromSelectionQuery(plotterParameters);

  unsigned int nodes = (DIMENSIONS == 3 ? _order  : 0 ) + 1;
  nodes *= (_order + 1) * (_order + 1);
  _tempSolution.resize(_solverUnknowns * nodes);
  _tempGradient.resize(DIMENSIONS * _solverUnknowns * nodes);
  assertion(_tempSolution.size()== _solverUnknowns * nodes);
  assertion(_tempGradient.size()==DIMENSIONS * _solverUnknowns * nodes);

  _resolution = 0;
  if (_plotterParameters.hasKey("resolution")) {
    _resolution = _plotterParameters.getValueAsIntOrDefault("resolution",0);
  }
  logInfo("init", "Plotting with resolution "<<_resolution);

  if(_slicer) {
    logInfo("init", "Plotting selection "<<_slicer->toString()<<" to Files "<<filename);
  }
}


void exahype::plotters::LimitingADERDG2CartesianVTKwithGradients::startPlotting( double time ) {
  _fileCounter++;

  assertion( _patchWriter==nullptr );

  if (_writtenUnknowns>0) {
    switch (_plotterType) {
    case PlotterType::BinaryVTK:
      _patchWriter =
          new tarch::plotter::griddata::blockstructured::PatchWriterUnstructured(
              new tarch::plotter::griddata::unstructured::vtk::VTKBinaryFileWriter());
      break;
    case PlotterType::ASCIIVTK:
      _patchWriter =
          new tarch::plotter::griddata::blockstructured::PatchWriterUnstructured(
              new tarch::plotter::griddata::unstructured::vtk::VTKTextFileWriter());
      break;
    case PlotterType::BinaryVTU:
      _patchWriter =
          new tarch::plotter::griddata::blockstructured::PatchWriterUnstructured(
              new tarch::plotter::griddata::unstructured::vtk::VTUBinaryFileWriter());
      break;
    case PlotterType::ASCIIVTU:
      _patchWriter =
          new tarch::plotter::griddata::blockstructured::PatchWriterUnstructured(
              new tarch::plotter::griddata::unstructured::vtk::VTUTextFileWriter());
      break;
    }

    _gridWriter                  = _patchWriter->createSinglePatchWriter();
    if ( true || _plotCells) {
      _cellDataWriter            = _patchWriter->createCellDataWriter("Q", _writtenUnknowns);
      _vertexDataWriter          = nullptr;

      _cellRefinementStatusWriter           = _patchWriter->createCellDataWriter("RefinementStatus", 1);
      _vertexRefinementStatusWriter         = nullptr;
      _cellPreviousRefinementStatusWriter   = _patchWriter->createCellDataWriter("PreviousRefinementStatus", 1);
      _vertexPreviousRefinementStatusWriter = nullptr;

      _timeStampVertexDataWriter   = nullptr;
    }
    else {
      _cellDataWriter            = _patchWriter->createCellDataWriter("Q", _writtenUnknowns);
      _vertexDataWriter          = _patchWriter->createVertexDataWriter("Q", _writtenUnknowns);

      _cellRefinementStatusWriter   = nullptr;
      _vertexRefinementStatusWriter = _patchWriter->createVertexDataWriter("RefinementStatus", 1);

      _cellPreviousRefinementStatusWriter   = nullptr;
      _vertexPreviousRefinementStatusWriter = _patchWriter->createVertexDataWriter("PreviousRefinementStatus", 1);

      _timeStampVertexDataWriter = _patchWriter->createVertexDataWriter("time", 1);
      assertion( _timeStampVertexDataWriter!=nullptr );
    }

    _timeStampCellDataWriter   = _patchWriter->createCellDataWriter("time", 1);
    assertion( _timeStampCellDataWriter!=nullptr );

    assertion( _patchWriter!=nullptr );
    assertion( _gridWriter!=nullptr );
  }

  _postProcessing->startPlotting( time );

  _time = time;
}


void exahype::plotters::LimitingADERDG2CartesianVTKwithGradients::finishPlotting() {
  _postProcessing->finishPlotting();

  if (_writtenUnknowns>0) {
    assertion( _patchWriter!=nullptr );
    assertion( _gridWriter!=nullptr );

    _gridWriter->close();
    if (_timeStampVertexDataWriter!=nullptr)           _timeStampVertexDataWriter->close();
    if (_timeStampCellDataWriter!=nullptr)             _timeStampCellDataWriter->close();
    if (_vertexDataWriter!=nullptr)                    _vertexDataWriter->close();
    if (_cellDataWriter!=nullptr)                      _cellDataWriter->close();
    if (_cellRefinementStatusWriter!=nullptr)           _cellRefinementStatusWriter->close();
    if (_vertexRefinementStatusWriter!=nullptr)         _vertexRefinementStatusWriter->close();
    if (_cellPreviousRefinementStatusWriter!=nullptr)   _cellPreviousRefinementStatusWriter->close();
    if (_vertexPreviousRefinementStatusWriter!=nullptr) _vertexPreviousRefinementStatusWriter->close();

    std::ostringstream snapshotFileName;
    snapshotFileName << _filename << "-" << _fileCounter;

    switch (_plotterType) {
      case PlotterType::BinaryVTK:
        break;
      case PlotterType::ASCIIVTK:
        break;
      case PlotterType::BinaryVTU:
        _timeSeriesWriter.addSnapshot( snapshotFileName.str(), _time);
        _timeSeriesWriter.writeFile(_filename);
        break;
      case PlotterType::ASCIIVTU:
        _timeSeriesWriter.addSnapshot( snapshotFileName.str(), _time);
        _timeSeriesWriter.writeFile(_filename);
        break;
    }

    const bool hasBeenSuccessful =
      _patchWriter->writeToFile(snapshotFileName.str());
    if (!hasBeenSuccessful) {
      exit(-1);
    }
  }

  if (_vertexDataWriter!=nullptr)                  delete _vertexDataWriter;
  if (_cellDataWriter!=nullptr)                    delete _cellDataWriter;
  if (_timeStampVertexDataWriter!=nullptr)         delete _timeStampVertexDataWriter;
  if (_timeStampCellDataWriter!=nullptr)           delete _timeStampCellDataWriter;
  if (_cellRefinementStatusWriter!=nullptr)           delete _cellRefinementStatusWriter;
  if (_vertexRefinementStatusWriter!=nullptr)         delete _vertexRefinementStatusWriter;
  if (_cellPreviousRefinementStatusWriter!=nullptr)   delete _cellPreviousRefinementStatusWriter;
  if (_vertexPreviousRefinementStatusWriter!=nullptr) delete _vertexPreviousRefinementStatusWriter;
  if (_gridWriter!=nullptr)                        delete _gridWriter;
  if (_patchWriter!=nullptr)                       delete _patchWriter;

  _vertexDataWriter                  = nullptr;
  _cellDataWriter                    = nullptr;
  _patchWriter                       = nullptr;
  _timeStampVertexDataWriter         = nullptr;
  _timeStampCellDataWriter           = nullptr;
  _cellRefinementStatusWriter           = nullptr;
  _vertexRefinementStatusWriter         = nullptr;
  _cellPreviousRefinementStatusWriter   = nullptr;
  _vertexPreviousRefinementStatusWriter = nullptr;
  _gridWriter                        = nullptr;
}



exahype::plotters::LimitingADERDG2CartesianVTKwithGradients::~LimitingADERDG2CartesianVTKwithGradients() {
}


void exahype::plotters::LimitingADERDG2CartesianVTKwithGradients::writeTimeStampDataToADERDGPatch( double timeStamp, int vertexOrCellIndex ) {
  if ( _plotCells && _writtenUnknowns > 0 ) {
    dfor(i,_order) {
      _timeStampCellDataWriter->plotCell(vertexOrCellIndex, timeStamp);
      vertexOrCellIndex++;
    }
  } else if ( _writtenUnknowns>0) {
    dfor(i,_order+1) {
      _timeStampVertexDataWriter->plotVertex(vertexOrCellIndex, timeStamp);
      vertexOrCellIndex++;
    }
  }
}


void exahype::plotters::LimitingADERDG2CartesianVTKwithGradients::plotVertexData(
  int firstVertexIndex,
  const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
  const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
  double* u,
  double timeStamp,
  const int RefinementStatusAsInt,
  const int previousRefinementStatusAsInt
) {
  assertion( _vertexDataWriter!=nullptr || _writtenUnknowns==0 );

  double* interpoland = new double[_solverUnknowns];
  double* value       = _writtenUnknowns==0 ? nullptr : new double[_writtenUnknowns];

  dfor(i,_order+1) {
    for (int unknown=0; unknown < _solverUnknowns; unknown++) {
      interpoland[unknown] = 0.0;
      dfor(ii,_order+1) { // Gauss-Legendre node indices
        int iGauss = peano::utils::dLinearisedWithoutLookup(ii,_order + 1);
        interpoland[unknown] += kernels::legendre::equidistantGridProjector[_order][ii(1)][i(1)] *
                 kernels::legendre::equidistantGridProjector[_order][ii(0)][i(0)] *
                 #if DIMENSIONS==3
                 kernels::legendre::equidistantGridProjector[_order][ii(2)][i(2)] *
                 #endif
                 u[iGauss * _solverUnknowns + unknown];
        assertion3(interpoland[unknown] == interpoland[unknown], offsetOfPatch, sizeOfPatch, iGauss);
      }
      if ( !std::isfinite(interpoland[unknown]) ) {
        logError("plotVertexData(...)","plotted value not finite.");
        std::abort();
      }
    }

    assertion(sizeOfPatch(0)==sizeOfPatch(1));
    _postProcessing->mapQuantities(
      offsetOfPatch,
      sizeOfPatch,
      offsetOfPatch + i.convertScalar<double>()* (sizeOfPatch(0)/(_order)),
      i,
      interpoland,
      value,
      timeStamp
    );

    if (_writtenUnknowns>0) {
      _vertexDataWriter->plotVertex(firstVertexIndex, value, _writtenUnknowns );
    }

    _vertexRefinementStatusWriter->plotVertex(firstVertexIndex, static_cast<double>(RefinementStatusAsInt));
    _vertexPreviousRefinementStatusWriter->plotVertex(firstVertexIndex, static_cast<double>(previousRefinementStatusAsInt));

    firstVertexIndex++;
  }

  if (interpoland!=nullptr)  delete[] interpoland;
  if (value!=nullptr)        delete[] value;
}


void exahype::plotters::LimitingADERDG2CartesianVTKwithGradients::plotCellData(
  int firstCellIndex,
  const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
  const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
  double* u,
  double timeStamp,
  const int RefinementStatusAsInt,
  const int previousRefinementStatusAsInt
) {
  assertion( _cellDataWriter!=nullptr || _writtenUnknowns==0 );

  double* interpoland = new double[_solverUnknowns];
  double* value       = _writtenUnknowns==0 ? nullptr : new double[_writtenUnknowns];

  dfor(i,_order) {
    for (int unknown=0; unknown < _solverUnknowns; unknown++) {
      interpoland[unknown] = kernels::legendre::interpolate(
        offsetOfPatch.data(),
        sizeOfPatch.data(),
        (offsetOfPatch + (i.convertScalar<double>()+0.5)* (sizeOfPatch(0)/(_order))).data(),
        _solverUnknowns,
        unknown,
        _order,
        u
      );
    }

    assertion(sizeOfPatch(0)==sizeOfPatch(1));
    _postProcessing->mapQuantities(
      offsetOfPatch,
      sizeOfPatch,
      offsetOfPatch + (i.convertScalar<double>()+0.5)* (sizeOfPatch(0)/(_order)),
      i,
      interpoland,
      value,
      timeStamp
    );

    if (_writtenUnknowns>0) {
      _cellDataWriter->plotCell(firstCellIndex, value, _writtenUnknowns );
    }

    _cellRefinementStatusWriter->plotCell(firstCellIndex, static_cast<double>(RefinementStatusAsInt));
    _cellPreviousRefinementStatusWriter->plotCell(firstCellIndex, static_cast<double>(previousRefinementStatusAsInt));

    firstCellIndex++;
  }

  if (interpoland!=nullptr)  delete[] interpoland;
  if (value!=nullptr)        delete[] value;
}

void exahype::plotters::LimitingADERDG2CartesianVTKwithGradients::plotPatch(const int solverNumber,solvers::Solver::CellInfo& cellInfo) {
  // look up ADER-DG solver
  solvers::ADERDGSolver*        aderdgSolver = nullptr;
  solvers::FiniteVolumesSolver* fvSolver = nullptr;
  switch ( solvers::RegisteredSolvers[solverNumber]->getType() ) {
  case solvers::Solver::Type::ADERDG:
    aderdgSolver = static_cast<solvers::ADERDGSolver*>( solvers::RegisteredSolvers[solverNumber] );
    break;
  case solvers::Solver::Type::LimitingADERDG:
    aderdgSolver = static_cast<solvers::LimitingADERDGSolver*>( solvers::RegisteredSolvers[solverNumber] )->getSolver().get();
    fvSolver     = static_cast<solvers::LimitingADERDGSolver*>( solvers::RegisteredSolvers[solverNumber] )->getLimiter().get();
    break;
  default:
    logError("plotPatch(...)","Encountered unexpected solver type: "<<solvers::Solver::toString(solvers::RegisteredSolvers[solverNumber]->getType()));
    std::abort();
    break;
  }

  const int element = cellInfo.indexOfADERDGCellDescription(solverNumber);
  auto& solverPatch  = cellInfo._ADERDGCellDescriptions[element];

  if ( solverPatch.getType()==exahype::solvers::ADERDGSolver::CellDescription::Type::Leaf ) {
    int refinementStatus         = solverPatch.getRefinementStatus();
    int previousRefinementStatus = solverPatch.getPreviousRefinementStatus();

    // ignore limiter status on coarser mesh levels
    assertion(static_cast<unsigned int>(solverPatch.getSolverNumber())<exahype::solvers::RegisteredSolvers.size());
    if (solverPatch.getLevel()<exahype::solvers::RegisteredSolvers[solverPatch.getSolverNumber()]->getMaximumAdaptiveMeshLevel()) {
      refinementStatus         = 0;
      previousRefinementStatus = 0;
    }

    if( refinementStatus < aderdgSolver->getMaxRefinementStatus()-1 ) {  // TODO(Dominic): Plot FVM solution instead if < Troubled-1
      double* solution = static_cast<double*>(solverPatch.getSolution());
      const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch = solverPatch.getOffset();

      const int subcellsPerDim = tarch::la::aPowI(_resolution,3);
      const tarch::la::Vector<DIMENSIONS, double>& subcellSize = solverPatch.getSize() / static_cast<double>(subcellsPerDim);

      dfor(subcellIndex,subcellsPerDim) {
        tarch::la::Vector<DIMENSIONS, double> subcellOffset = offsetOfPatch;
        double* u = solution;
        if ( subcellsPerDim > 1 ) {
          u = _tempSolution.data();
          for (int d=0; d<DIMENSIONS; d++) {
            subcellOffset[d] = offsetOfPatch[d] + subcellSize[d] * subcellIndex[d];
          }
          aderdgSolver->volumeUnknownsProlongation(u,solution,0,_resolution,subcellIndex);
        }

        plotADERDGPatch(
            subcellOffset,
            subcellSize,
            u,
            solverPatch.getTimeStamp(),
            refinementStatus,
            previousRefinementStatus);
      }
    } else {
      auto& limiterPatch = cellInfo._FiniteVolumesCellDescriptions[solverNumber];
      plotFiniteVolumesPatch
          (solverPatch.getOffset(),solverPatch.getSize(),
          static_cast<double*>(limiterPatch.getSolution()),
          solverPatch.getTimeStamp(),
          fvSolver->getNodesPerCoordinateAxis(),
          refinementStatus,previousRefinementStatus);
    }
  }
}


void exahype::plotters::LimitingADERDG2CartesianVTKwithGradients::plotADERDGPatch(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    double* u,
    double timeStamp,
    const int RefinementStatusAsInt,
    const int previousRefinementStatusAsInt) {
  if (!_slicer || _slicer->isPatchActive(offsetOfPatch, sizeOfPatch)) {
    assertion( _writtenUnknowns==0 || _patchWriter!=nullptr );
    assertion( _writtenUnknowns==0 || _gridWriter!=nullptr );

    std::pair<int,int> vertexAndCellIndex(0,0);
    if (_writtenUnknowns>0) {
      vertexAndCellIndex = _gridWriter->plotPatch(offsetOfPatch, sizeOfPatch, _order);
    }

    if (_plotCells) {
      writeTimeStampDataToADERDGPatch( timeStamp, vertexAndCellIndex.second );

      assertion( _writtenUnknowns==0 || _timeStampCellDataWriter!=nullptr );
      plotCellData( vertexAndCellIndex.second, offsetOfPatch, sizeOfPatch, u, timeStamp, RefinementStatusAsInt, previousRefinementStatusAsInt );
    }
    else {
      writeTimeStampDataToADERDGPatch( timeStamp, vertexAndCellIndex.first );

      assertion( _writtenUnknowns==0 || _timeStampVertexDataWriter!=nullptr );
      plotVertexData( vertexAndCellIndex.first, offsetOfPatch, sizeOfPatch, u, timeStamp, RefinementStatusAsInt, previousRefinementStatusAsInt );
    }
  }
}

void exahype::plotters::LimitingADERDG2CartesianVTKwithGradients::plotFiniteVolumesPatch(
  const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
  const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
  double*                                      u,
  double                                       timeStamp,
  const int                                    numberOfCellsPerAxis,
  int                                          RefinementStatusAsInt,
  int                                          previousRefinementStatusAsInt) {
  if (!_slicer || _slicer->isPatchActive(offsetOfPatch, sizeOfPatch)) {
    logDebug("plotPatch(...)","offset of patch: "<<offsetOfPatch
    <<", size of patch: "<<sizeOfPatch
    <<", time stamp: "<<timeStamp);

    assertion( _writtenUnknowns==0 || _patchWriter!=nullptr );
    assertion( _writtenUnknowns==0 || _gridWriter!=nullptr );
    assertion( _writtenUnknowns==0 || _timeStampCellDataWriter!=nullptr );

    int cellIndex = _writtenUnknowns==0 ? -1 : _gridWriter->plotPatch(offsetOfPatch, sizeOfPatch, numberOfCellsPerAxis).second;

    double* sourceValue = new double[_solverUnknowns];
    double* value       = _writtenUnknowns==0 ? nullptr : new double[_writtenUnknowns];

    dfor(i,numberOfCellsPerAxis+_ghostLayerWidth) {
      if (tarch::la::allSmaller(i,numberOfCellsPerAxis+_ghostLayerWidth)
          && tarch::la::allGreater(i,_ghostLayerWidth-1)) {
        if (_writtenUnknowns>0) {
          _timeStampCellDataWriter->plotCell(cellIndex, timeStamp);
        }

        for (int unknown=0; unknown < _solverUnknowns; unknown++) {
          sourceValue[unknown] =
            u[peano::utils::dLinearisedWithoutLookup(i,numberOfCellsPerAxis+2*_ghostLayerWidth)*_solverUnknowns+unknown];
        } // !!! Be aware of the "2*_ghostLayerWidth" !!!

        assertion(sizeOfPatch(0)==sizeOfPatch(1));

        _postProcessing->mapQuantities(
          offsetOfPatch,
          sizeOfPatch,
          offsetOfPatch + (i-_ghostLayerWidth).convertScalar<double>()* (sizeOfPatch(0)/(numberOfCellsPerAxis)),
          i-_ghostLayerWidth,
          sourceValue,
          value,
          timeStamp
        );

        if (_writtenUnknowns>0) {
          _cellDataWriter->plotCell(cellIndex, value, _writtenUnknowns);
        }

        _cellRefinementStatusWriter->plotCell(cellIndex, static_cast<double>(RefinementStatusAsInt));
        _cellPreviousRefinementStatusWriter->plotCell(cellIndex, static_cast<double>(previousRefinementStatusAsInt));

        cellIndex++;
      }
    }

    delete[] sourceValue;
    delete[] value;
  }
}
