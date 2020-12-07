#pragma once
#include <string>
#include <vector>

class STATE
{
public:
	double u;
	double s;
	double ss;
};

class CELL
{
public:
	// Data Members
	double xL;
	double xR;
	double uM;
	double uL;
	double uR;
	double sL;
	double sR;
	double ssL;
	double ssR;
};

class MESH
{
public:
	// Data Members
	int p;
	int NumCells;
	double Dx;
	CELL *cells;
};