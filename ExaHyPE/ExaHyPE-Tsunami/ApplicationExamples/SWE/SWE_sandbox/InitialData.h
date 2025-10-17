#ifndef __InitialData_CLASS_HEADER__
#define __InitialData_CLASS_HEADER__

namespace easi {
  class YAMLParser;
  class Component;
};
class AsagiReader;

class InitialData {
 public:
  InitialData();
  InitialData(int scenario, char* filename);  
  ~InitialData();
  void getInitialData(const double* const x,double* Q);
  
 private:
  int scenario;
  easi::YAMLParser* parser;
  easi::Component* model;
  AsagiReader* asagiReader;

  void readAsagiData(const double* const x,double* Q);
  void readAsagiData_nobath(const double* const x,double* Q);
};

#endif // __InitialData_CLASS_HEADER__
