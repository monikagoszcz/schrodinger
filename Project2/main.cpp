/*
 * main.cpp
 *
 *  Created on: Oct 5, 2016
 *      Author: mongos
 */

#include "simulation.h"

int main(int argc, char*  argv[])
{
	std::string number = argv[1];
	Parameters parameters = getParameters("data"+ number + ".txt");
	
	State state = {};

	state.H.resize(parameters.N + 1);
	state.ro.resize(parameters.N + 1);
	state.fi.resize(parameters.N + 1);
	setInitialState(parameters, state);

	std::ofstream outputFile;
    outputFile.open("output" + number +".txt");

	std::ofstream outputResonanceFile;
	outputResonanceFile.open("outputR" + number + ".txt");

    simulate(parameters, state, outputFile, outputResonanceFile);

    return 0;
}


