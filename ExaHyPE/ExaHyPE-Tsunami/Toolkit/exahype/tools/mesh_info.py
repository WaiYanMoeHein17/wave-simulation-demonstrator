#author: Dominic E. Charrier

import os,sys,math

def getCoarsestMaximumMeshSizeOfAllSolvers(solverSpec):
  """
  Find the coarsest maximum mesh size the user has specified
  for a solver.
  """
  result = 1e20
  for solver in solverSpec:
    result = min( result, solver["maximum_mesh_size"] )
  return result

class MeshInfo:
  def __init__(self,domainOffset,domainSize,userMeshSize,outsideCellsLeft,outsideCellsRight,ranksPerDimension):
    """
    Constructor. All arguments are required/

    :type vector: list[float]
    :param vector domainOffset:       offset of the domain   
    :param vector domainSize:         size of the domain
    :param float userMeshSize:        the mesh size the user at least wants
    :param float outsideCellsLeft:    number of cells which should be placed outside to the left of the domain
    :param float outsideCellsRight:   number of cells which should be placed outside to the right of the domain
    :param bool  ranksPerDimension:   If the user wants to use a specific number of cells outside of the domain. 
                                      If this is a multiple of 3, the outside cells are chosen as 2.
                                      In this case, one outside cell is placed left and the other right to the domain per dimension.
                                      Overrules the other scaling options.
    """
    self.domainOffset       = domainOffset
    self.domainSize         = domainSize
    self.unscaledDomainSize = domainSize
    self.userMeshSize       = userMeshSize
    self.dim                = len(self.domainOffset)

    self.boundingBoxOutsideCells      = outsideCellsLeft + outsideCellsRight # is overwritten if ranksPerDimension is true
    self.boundingBoxOutsideCellsLeft  = outsideCellsLeft  # is overwritten if ranksPerDimension is true
    self.boundingBoxOutsideCellsRight = outsideCellsRight # is overwritten if ranksPerDimension is true
    self.ranksPerDimension            = ranksPerDimension 

    # deduced quantities
    self.scaleBoundingBox       = (outsideCellsRight+outsideCellsLeft)>0 or ranksPerDimension>0
    self.boundingBoxMeshSize    = -1
    self.boundingBoxMeshLevel   = -1
    self.boundingBoxMeshCells   = -1
    self.boundingBoxSize        = -1
    self.boundingBoxOffset      = -1
    self.boundingBoxInsideCells = -1

  def computeCoarsestMeshSizeAndLevel(self):
    """
    Enlarge non-cubic domains such that they are exactly
    resolved by a grid with mesh size meshSize.
    """
    largestExtent = max(self.domainSize)

    level = 0 
    currentMeshSize = 1e20
    while currentMeshSize>self.userMeshSize:
      currentMeshSize = float(largestExtent) / 3.0**level
      level += 1
    level -= 1 # currentMeshSize was computed with previous level (= new level -1 )
    return currentMeshSize, level

  def determineScaledDomainSize(self,meshSize):
    """
    Enlarge non-cubic domains such that they are exactly
    resolved by a grid with mesh size meshSize.
    """
    scaledDomainSize=[0 for i in range(0,self.dim)]    
    for d in range(0,self.dim):
      scaledDomainSize[d] = math.ceil(self.domainSize[d] / meshSize - 1e-11) * meshSize;
    return scaledDomainSize

  def deduceCoarseGridInformation(self):
    """
    From the solver specification and the bounding box scaling options,
    deduce how the coarse grid will look like.
    """
    self.boundingBoxSize                = max(self.domainSize);
    unscaledMeshSize, unscaledMeshLevel = self.computeCoarsestMeshSizeAndLevel();    
    self.boundingBoxMeshLevel           = unscaledMeshLevel;
    self.boundingBoxOffset              = self.domainOffset.copy()

    # scale bounding box
    if self.scaleBoundingBox:
      levelLB = 0
      powerOfThreeRanks = True
      if self.ranksPerDimension>0:
        levelLB = int(math.ceil( math.log(self.ranksPerDimension)/math.log(3) - 1e-9 ))
        powerOfThreeRanks = self.ranksPerDimension is 3**levelLB

      maxDomainExtent = max(self.domainSize)
      self.boundingBoxMeshSize = -1
      boundingBoxScaling  = 0
      level = unscaledMeshLevel; # level=0 means a single cell

      self.boundingBoxInsideCells = 0
      while level <= levelLB+1 or\
            self.boundingBoxOutsideCellsRight < 0 or\
            self.boundingBoxMeshSize < 0 or\
            self.boundingBoxMeshSize > self.userMeshSize:
        self.boundingBoxMeshCells = 3**level;

        self.boundingBoxOutsideCells = self.boundingBoxOutsideCellsLeft + self.boundingBoxOutsideCellsRight
        if not powerOfThreeRanks: # overwrite
          self.boundingBoxOutsideCells = int(round(self.boundingBoxMeshCells * (1.0-self.ranksPerDimension/3.0**levelLB))) # must be >= left

        self.boundingBoxInsideCells = int(self.boundingBoxMeshCells -self. boundingBoxOutsideCells)
        boundingBoxScaling          = float(self.boundingBoxMeshCells) / float(self.boundingBoxInsideCells)
        self.boundingBoxSize        = boundingBoxScaling * maxDomainExtent
        self.boundingBoxMeshSize    = self.boundingBoxSize / self.boundingBoxMeshCells
        self.boundingBoxOutsideCellsRight = self.boundingBoxOutsideCells - self.boundingBoxOutsideCellsLeft
        level += 1
     
      self.boundingBoxMeshLevel = level - 1; # decrement result since bounding box was computed using level-1
      for d in range(0,self.dim):
        self.boundingBoxOffset[d] -= self.boundingBoxOutsideCellsLeft*self.boundingBoxMeshSize
    
    else: # not scale bounding box
      self.boundingBoxMeshSize    = unscaledMeshSize
      self.boundingBoxMeshCells   = 3**self.boundingBoxMeshLevel;
      self.boundingBoxInsideCells = self.boundingBoxMeshCells

    self.domainSize = self.determineScaledDomainSize(self.boundingBoxMeshSize)
