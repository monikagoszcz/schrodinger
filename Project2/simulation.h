/*
 * argon.h
 *
 *  Created on: Oct 5, 2016
 *      Author: mongos
 */

#ifndef _simulation_h_
#define _simulation_h_

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

struct Complex
{
	double re;
	double im;
};

struct ResonantParameters
{
	std::vector<double> omega;
	std::vector<double> maxEnergy;
};

struct Parameters
{
    int N;
	double kappa;
	double omega;
	double dtau;
	double dx;
	int Sout;
	int Stot;
	double n; // if = 1 =, then it is the ground state
};

struct State
{
	double N; //norma
    double x; //srednie polozenie
    double E; //energia
	double t;
	std::vector<double> ro;
	std::vector<Complex> H;
	std::vector<Complex> fi;
};

Parameters getParameters(std::string inputFileName);

void setInitialState(const Parameters &Parameters, State & state);
void setInitialFi(const Parameters &parameters, std::vector<Complex> & fi);
void setHamiltonian(const Parameters &parameters, State & state, double &t);

void simulate(const Parameters &Parameters, State &state, std::ofstream &outputFileState);

void updateState(const Parameters &Parameters, State & state);
void setStateParameters(const Parameters &Parameters, State & state);

void outputState(State &state, std::ofstream &outputFile);
void outputMaxEnergy(ResonantParameters &res, std::ofstream &outputFileResonance);

void findResonanceCurve(Parameters &Parameters, State &state, std::ofstream &outputFileResonance);
void simulate(const Parameters &Parameters, State &state, std::vector<double> & Energies, std::ofstream & outputFileChar);
double maxValue(std::vector<double> & vect);


#endif 

