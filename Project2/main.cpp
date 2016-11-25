/*
 * main.cpp
 *
 *  Created on: Oct 5, 2016
 *      Author: mongos
 */

#include "simulation.h"
#include "coor.h"

int main(int /* argc */, char* /* argv */[])
{
	Parameters parameters = getParameters("data.txt");
	
	State state = {};

	state.H.resize(parameters.N + 1);
	state.ro.resize(parameters.N + 1);
	state.fi.resize(parameters.N + 1);
	setInitialState(parameters, state);

	std::ofstream outputFile;
    outputFile.open("output.txt");

    simulate(parameters, state, outputFile);

	std::ofstream outputResonanceFile;
	outputResonanceFile.open("outputR.txt");

	findResonanceCurve(parameters, state, outputResonanceFile);


    return 0;
}


