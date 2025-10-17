#include "iodir.h"

std::string get_input(){
    const char* env_path = std::getenv("SHARED_DIR");
    if(env_path == NULL){
        env_path = "/tmp/";
        std::cout << "Shared directory path not set" <<std::endl;
    }
    std::string input(env_path);
    return input+"inputs.txt";
}

std::string get_output(std::string probe_name){
    const char* env_path = std::getenv("SHARED_DIR");
    if(env_path == NULL){
        env_path = "/tmp/";
        std::cout << "Shared directory path not set" <<std::endl;
    }
    std::string output(env_path);
    return output+probe_name+"outputs.txt";
}
