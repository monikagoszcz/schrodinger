/*
 * argon.cpp
 *
 *  Created on: Oct 5, 2016
 *      Author: mongos
 */

#include "simulation.h"
#include "Coor.h"

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstring>

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
		state.H[k].re = {-0.5 * (state.fi[k+1].re + state.fi[k].re - 2 * state.fi[k-1].re) / (parameters.dx  * parameters.dx)
						+ parameters.kappa * (xk - .5) * state.fi[k].re * sin(parameters.omega * t)};
		state.H[k].im = { -0.5 * (state.fi[k + 1].im + state.fi[k].im - 2 * state.fi[k - 1].im) / (parameters.dx  * parameters.dx)
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
		state.fi[k].im = state.fi[k].im + .5 * state.H[k].re * parameters.dtau;
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

