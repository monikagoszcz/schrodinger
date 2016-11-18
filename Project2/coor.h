#ifndef _coor_h_
#define _coor_h_

#include <iostream>

struct Coor
{
	double x;
	double y;
	double z;
};

double calcVectorModulus(Coor vector);
Coor subtractVectors(Coor v1, Coor v2);
Coor addVectors(Coor v1, Coor v2);
Coor calcOppositeVector(Coor v);
#endif 

