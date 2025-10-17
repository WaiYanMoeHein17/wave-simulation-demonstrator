#ifndef ASAGIREADER_CLASS_HEADER
#define ASAGIREADER_CLASS_HEADER

#include <asagi.h>
#include <easi/util/AsagiReader.h>
#include <iostream>

#ifdef USE_ASAGI

enum NUMACache_Mode
{
	NUMA_OFF, NUMA_ON, NUMA_CACHE
};

class AsagiReader : public easi::AsagiReader{
 private:
  const std::string m_envPrefix;
 public:
  AsagiReader(const char* envPrefix) : m_envPrefix(envPrefix) {}
  virtual asagi::Grid* open(char const* file, char const* varname);
  virtual unsigned numberOfThreads() const { return 1; }
};

asagi::Grid* AsagiReader::open(char const* file, char const* varname){

  asagi::Grid* grid = asagi::Grid::create();
  
  //  grid->setParam(numberOfThreads());
  grid->setParam("NUMA_COMMUNICATION", "OFF",0);
  
  grid->setParam("GRID", "CACHE",0);

  grid->setParam("VALUE_POSITION","VERTEX_CENTERED",0);
  grid->setParam("VARIABLE",varname,0);

  grid->setParam("BLOCK_SIZE_0", "64",0);
  grid->setParam("BLOCK_SIZE_1", "64",0);
  grid->setParam("BLOCK_SIZE_2", "64",0);
  grid->setParam("CACHE_SIZE", "128",0);
  // grid->setParam("LEVEL", "0");

  asagi::Grid::Error err = grid->open(file,0);

  if(err != asagi::Grid::SUCCESS){
    std::cout << "Could not open "<< file << std::endl;
  }

  return grid;
}


#endif
#endif
