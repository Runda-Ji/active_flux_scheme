#pragma once
#include <string>
#include <vector>

class STATE
{
public:
	int p;
	int NumCells;
	double *u;
};

class CELL
{
public:
	// Data Members
	double xL;
	double xR;
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

class QUADRATURE
{
public:
	// Data Members
	int n;
	double *x;
	double *w;
};