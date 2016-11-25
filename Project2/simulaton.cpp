/*
 * argon.cpp
 *
 *  Created on: Oct 5, 2016
 *      Author: mongos
 */

#include "simulation.h"

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <algorithm>

using namespace std;

#ifndef _WIN32 /*Enables compilation not only on Windows*/
int fopen_s(FILE **f, const char *name, const char *mode) {
    int ret = 0;
    assert(f);
    *f = fopen(name, mode);
    if (!*f)
        ret = errno;
    return ret;
}
#define fscanf_s fscanf
#endif

/* constants */
namespace {
	const double PI = 3.141592653589793;
}

Parameters getParameters(std::string inputFileName)
{
    Parameters param = {};
	std::ifstream inputFile;
	inputFile.open(inputFileName);
 
		inputFile  >> param.N
		           >> param.kappa
                   >> param.omega
		           >> param.dtau
		           >> param.Sout
		           >> param.Stot
				   >> param.n;
		param.dx = 1 / static_cast<double>(param.N);
		param.omega = param.omega * PI * PI;

    return param;
}

void setInitialState(const Parameters &parameters, State &state)
{
    setInitialFi(parameters, state.fi);
    setHamiltonian(parameters, state, state.t);
}

void setInitialFi(const Parameters &parameters, vector<Complex> & fi)
{

    for(int k = 0; k < (parameters.N + 1); k++)
    {
		double xk = k * parameters.dx;
		fi[k].re = {sqrt(2)*sin(parameters.n*PI*xk)};
		fi[k].im = { 0 };
    }
}

void setHamiltonian(const Parameters &parameters, State & state, double &t)
{

	state.H[0].re = 0;
	state.H[0].im = 0;
	state.H[parameters.N].re = 0;
	state.H[parameters.N].im = 0;
	double xk = 0;
	for (int k = 1; k < parameters.N; k++)
	{
		xk += parameters.dx;
		state.H[k].re = {-0.5 * (state.fi[k+1].re + state.fi[k - 1].re - 2 * state.fi[k].re) / (parameters.dx  * parameters.dx)
						+ parameters.kappa * (xk - .5) * state.fi[k].re * sin(parameters.omega * t)};
		state.H[k].im = { -0.5 * (state.fi[k + 1].im + state.fi[k - 1].im - 2 * state.fi[k].im) / (parameters.dx  * parameters.dx)
			            + parameters.kappa * (xk - .5) * state.fi[k].im * sin(parameters.omega * t) };
	}
}


void simulate(const Parameters &parameters, State &state, ofstream & outputFileChar)
{
	state.t = 0;
    for(int s = 1; s <= (parameters.Stot); s++)
    {
        updateState(parameters, state);

		if (!(s % parameters.Sout))
		{
			setStateParameters(parameters, state);
			outputState(state, outputFileChar);
		}
    }
}

void updateState(const Parameters &parameters, State &state)
{
	for (int k = 0; k < (parameters.N + 1); k++)
    {
		state.fi[k].re = state.fi[k].re + .5 * state.H[k].im * parameters.dtau;
    }
	double t = (state.t + (.5 * parameters.dtau));
	setHamiltonian(parameters, state, t );

	for (int k = 0; k < (parameters.N + 1); k++)
	{
		state.fi[k].im = state.fi[k].im  -  state.H[k].re * parameters.dtau;
	}
	state.t += parameters.dtau;
	setHamiltonian(parameters, state, state.t);

	for (int k = 0; k < (parameters.N + 1); k++)
	{
		state.fi[k].re = state.fi[k].re + .5 * state.H[k].im * parameters.dtau;
	}
}

void setStateParameters(const Parameters &parameters, State &state)
{
	state.N = 0;
	state.x = 0;
	state.E = 0;
	double xk = 0; 
	for (int k = 0; k < (parameters.N + 1); k++)
	{
		state.ro[k] = state.fi[k].re * state.fi[k].re + state.fi[k].im * state.fi[k].im;
		state.N += state.ro[k];
		state.x += xk * state.ro[k];
		state.E += state.fi[k].re * state.H[k].re + state.fi[k].im * state.H[k].im;
		xk += parameters.dx;
	}
	state.N *= parameters.dx;
	state.x *= parameters.dx;
	state.E *= parameters.dx;
}


void outputState(State & state, ofstream &outputFile)
{
    outputFile << state.t << "\t" << state.N << "\t" << state.x << "\t" << state.E << "\t" << endl;
}

void outputMaxEnergy(ResonantParameters &res, ofstream &outputFileResonance)
{
	for (int i = 0; i < res.omega.size() - 1; i++)
	{
		outputFileResonance << res.omega[i] << "\t" << res.maxEnergy[i] << endl;
	}

}

void findResonanceCurve(Parameters &Parameters, State &state, ofstream &outputFileResonance) 
{
	ResonantParameters  res = {};
	int N = 30;
	res.omega.resize(N+1);
	res.maxEnergy.resize(N);
	double omegaBegin = 0.5 * Parameters.omega;
	double omegaEnd = 1.25 * Parameters.omega;
	double step = (omegaEnd - omegaBegin) / static_cast<double>(N);

	res.omega[0] = omegaBegin;

	for (int i = 0; i < N; i++)
	{
		string name = "output" + to_string(i) + ".txt";
		ofstream outputFile;
		outputFile.open(name);
		vector<double> Energies = {};
		Parameters.omega = res.omega[i];
		simulate(Parameters, state, Energies, outputFile);
		res.maxEnergy[i] = maxValue(Energies);
		res.omega[i + 1] = res.omega[i] + step;
	}
	outputMaxEnergy(res, outputFileResonance);
	
}

void simulate(const Parameters &Parameters, State &state, vector<double> & Energies, ofstream & outputFileChar)
{

	state.t = 0;
	for (int s = 1; s <= (Parameters.Stot); s++)
	{
		updateState(Parameters, state);

		if (!(s % Parameters.Sout))
		{
			setStateParameters(Parameters, state);
			Energies.push_back(state.E);
			outputState(state, outputFileChar);
		}
	}
}

double maxValue(vector<double> & vect)
{
	double max = vect[0];
	for (int i = 1; i < vect.size(); i++)
	{
		if (max < vect[i])
			max = vect[i];
	}
	return max;
}

