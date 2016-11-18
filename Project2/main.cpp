/*
 * main.cpp
 *
 *  Created on: Oct 5, 2016
 *      Author: mongos
 */

#include "simulation.h"
#include "coor.h"
#include "main.h"

int main(int /* argc */, char* /* argv */[])
{

	State state = {};
    Parameters parameters = getParameters("data.txt");

	setInitialState(parameters, state);

	std::ofstream outputFile;
    outputFile.open("output.txt");

    simulate(parameters, state, outputFile);

    return 0;
}


