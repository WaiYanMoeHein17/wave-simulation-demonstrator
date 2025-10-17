#include "InitialData.h"

#include "easi/YAMLParser.h"
#include "easi/ResultAdapter.h"
#include "reader/asagi_reader.h"
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include "iodir.h"

using namespace std;
///// 2D /////
std::vector<double> param = {0.0,0.0};
#define Dim2
#ifdef Dim2

InitialData::InitialData()
	: scenario(){
		std::cout << "Initialising with ASAGI" << std::endl;

        auto inputs = get_input();
        std::ifstream inputsfile(inputs);
        for (int i = 0; i < 2; i++) {
			inputsfile >> param[i];
		}
		inputsfile.close();
		std::cout << "Read inputs in exahype:" << param[0] << " and " << param[1] << std::endl;

		asagiReader = new AsagiReader("");
		parser = new easi::YAMLParser(3, asagiReader);
		model  = parser->parse("data.yaml");
	}

InitialData::InitialData(int a_scenario, char* filename)
	: scenario(a_scenario){
		std::cout << "Initialising with ASAGI" << std::endl;

        auto inputs = get_input();
        std::ifstream inputsfile(inputs);
        for (int i = 0; i < 2; i++) {
			inputsfile >> param[i];
		}
		inputsfile.close();
		std::cout << "Read inputs in exahype:" << param[0] << " and " << param[1] << std::endl;

		asagiReader = new AsagiReader("");
		parser = new easi::YAMLParser(3, asagiReader);
		model  = parser->parse(filename);
	}

InitialData::~InitialData(){
	delete asagiReader;
	delete parser;
	delete model;
}


void InitialData::readAsagiData(const double* const x,double* Q){
	double terrain[1];
	terrain[0] = 0.0;

	double water[1];
	water[0] = 0.0;

	easi::ArraysAdapter<double> adapter;
	adapter.addBindingPoint("b",terrain);

	easi::Query query(1,3);
	query.x(0,0)=x[0];
	query.x(0,1)=x[1];
	query.x(0,2)=0;
	model->evaluate(query,adapter);

	easi::ArraysAdapter<double> waterAdapter;
	waterAdapter.addBindingPoint("w",water);

	easi::Query query2(1,3);
	query2.x(0,0)=x[0];
	query2.x(0,1)=x[1];
	query2.x(0,2)=0;
	model->evaluate(query2,waterAdapter);

	double eta = 1.5; // water height

	// Water Height
	double sigma =  60.0; // width of initial wave
	double x_center = 400;
	double y_center = 300;
	// Q[0] = std::max(eta + 3.0 * exp(-((x[0] - x_center)*(x[0] - x_center) + (x[1] - y_center)*(x[1] - y_center)) / (2 * sigma*sigma)) - terrain[0], 0.0);
	Q[0] = std::max(water[0] - terrain[0], 0.0);
	// Q[0] = std::max(water[0] - terrain[0], 0.0);
	// Velocity
	Q[1] = 0.0;
	Q[2] = 0.0;
	// Bathymetry
	// Q[3] = terrain[0] > waterHeight ? 0.0 : -terrain[0];
	// Q[3] = terrain[0];
	Q[3]= water[0];
}

void InitialData::getInitialData(const double* const x,double* Q) {
	readAsagiData(x, Q);
}

#endif
