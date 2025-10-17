#include <iostream>

#include <string>

#include "umbridge.h"

#include <chrono>
#include <thread>
#include <iomanip>
#include <stdlib.h>

class TsunamiModel : public umbridge::Model {
public:

  TsunamiModel(int ranks)
   : Model("forward"), ranks(ranks)
  {
    char const* shared_dir_cstr = std::getenv("SHARED_DIR");
    if ( shared_dir_cstr == NULL ) {
      std::cerr << "Environment variable SHARED_DIR not set!" << std::endl;
      exit(-1);
    }
    shared_dir = std::string(shared_dir_cstr);
  }

  std::vector<std::size_t> GetInputSizes(const json& config) const override {
    return { 2 };
  }

  std::vector<std::size_t> GetOutputSizes(const json& config) const override {
    return {4};
  }

  std::vector<std::vector<double>> Evaluate(std::vector<std::vector<double>> const& inputs, json config) override {
    // TODO add any configurations you need like this
    //std::cout << "Reading options" << std::endl;
    //bool vtk_output = config.value("vtk_output", true);

    std::ofstream inputsfile (shared_dir + "inputs.txt");
    typedef std::numeric_limits<double> dl;
    inputsfile << std::fixed << std::setprecision(dl::digits10);
    for (int i = 0; i < inputs[0].size(); i++) {
      inputsfile << inputs[0][i] << std::endl;
    }
    inputsfile.close();

    int status;
    std::string cmd = "cd /ExaHyPE-Tsunami/ApplicationExamples/SWE/SWE_sandbox && mpirun --allow-run-as-root -x LD_LIBRARY_PATH -x SHARED_DIR -n " + std::to_string(ranks) + " ./ExaHyPE-SWE ../SWE_sandbox.exahype2 && cp output/* /output/";
    status = system(cmd.c_str());
    std::cout << "Exahype exit status " << status << std::endl;

    std::vector<std::vector<double>> output = {{0,0,0,0}};
    return output; // Don't need an output
  }

  bool SupportsEvaluate() override {
    return true;
  }

private:
  int ranks;
  std::string shared_dir;
};

int main(){

  char const* port_cstr = std::getenv("PORT");
  if ( port_cstr == NULL ) {
    std::cerr << "Environment variable PORT not set!" << std::endl;
    exit(-1);
  }
  const int port = atoi(port_cstr);

  char const* ranks_cstr =  std::getenv("RANKS");
  if ( ranks_cstr == NULL ) {
    std::cerr << "Environment variable RANKS not set!" << std::endl;
    exit(-1);
  }
  const int ranks = atoi(ranks_cstr);

  std::cout << "Running on number of ranks: " << ranks << std::endl;

  TsunamiModel model(ranks);
  std::vector<umbridge::Model*> models {&model};
  umbridge::serveModels(models, "0.0.0.0", port);

  return 0;
}
