void FindBasis(int p, double xi, double *phi)
{
	double xi2 = xi*xi;
	double xi3 = xi*xi*xi;
	double xi4 = xi*xi*xi*xi;
	double xi5 = xi*xi*xi*xi*xi;
	double xi6 = xi*xi*xi*xi*xi*xi;
	switch (p)
	{
	case 3: // AF3: uL uM uR
		phi[0] = 3*xi2 - 4*xi + 1;
		phi[1] = -6*xi2 + 6*xi;
		phi[2] = 3*xi2 - 2*xi;
		break;
	case 5: // AF5: uL sL uM uR sR
		phi[0] = -15*xi4 + 32*xi3 - 18*xi2 + 1;
		phi[1] = -5*xi4/2 + 6*xi3 - 9*xi2/2 + xi;
		phi[2] = 30*xi4 - 60*xi3 + 30*xi2;
		phi[3] = -15*xi4 + 28*xi3 - 12*xi2;
		phi[4] = 5*xi4/2 - 4*xi3 + 3*xi2/2;
		break;
	case 7: // AF7: uL sL ssL uM uR sR ssR
		phi[0] = 70*xi6 - 216*xi5 + 225*xi4 - 80*xi3 + 1;
		phi[1] = 14*xi6 - 45*xi5 + 50*xi4 - 20*xi3 + xi;
		phi[2] = 7*xi6/6 - 4*xi5 + 5*xi4 - 8*xi3/3 + xi2/2;
		phi[3] = -140*xi6 + 420*xi5 - 420*xi4 + 140*xi3;
		phi[4] = 70*xi6 - 204*xi5 + 195*xi4 - 60*xi3;
		phi[5] = -14*xi6 + 39*xi5 - 35*xi4 + 10*xi3;
		phi[6] = 7*xi6/6 - 3*xi5 + 5*xi4/2 - 2*xi3/3;
		break;
	}
}
void FindDiffBasis(int p, double xi, double *phi)
{
	double xi2 = xi*xi;
	double xi3 = xi*xi*xi;
	double xi4 = xi*xi*xi*xi;
	double xi5 = xi*xi*xi*xi*xi;
	switch (p)
	{
	case 3: // AF3: uL uM uR
		phi[0] = 6*xi - 4;
		phi[1] = 6 - 12*xi;
		phi[2] = 6*xi - 2;
		break;
	case 5: // AF5: uL sL uM uR sR
		phi[0] = -60*xi3 + 96*xi2 - 36*xi;
		phi[1] = -10*xi3 + 18*xi2 - 9*xi + 1;
		phi[2] = 120*xi3 - 180*xi2 + 60*xi;
		phi[3] = -60*xi3 + 84*xi2 - 24*xi;
		phi[4] = 10*xi3 - 12*xi2 + 3*xi;
		break;
	case 7: // AF7: uL sL ssL uM uR sR ssR
		phi[0] = 420*xi5 - 1080*xi4 + 900*xi3 - 240*xi2;
		phi[1] = 84*xi5 - 225*xi4 + 200*xi3 - 60*xi2 + 1;
		phi[2] = 7*xi5 - 20*xi4 + 20*xi3 - 8*xi2 + xi;
		phi[3] = -840*xi5 + 2100*xi4 - 1680*xi3 + 420*xi2;
		phi[4] = 420*xi5 - 1020*xi4 + 780*xi3 - 180*xi2;
		phi[5] = -84*xi5 + 195*xi4 - 140*xi3 + 30*xi2;
		phi[6] = 7*xi5 - 15*xi4 + 10*xi3 - 2*xi2;
		break;
	}
}
void FindDiff2Basis(int p, double xi, double *phi)
{
	double xi2 = xi*xi;
	double xi3 = xi*xi*xi;
	double xi4 = xi*xi*xi*xi;
	switch (p)
	{
	case 3: // AF3: uL uM uR
		phi[0] = 6;
		phi[1] = -12;
		phi[2] = 6;
		break;
	case 5: // AF5: uL sL uM uR sR
		phi[0] = -180*xi2 + 192*xi - 36;
		phi[1] = -30*xi2 + 36*xi - 9;
		phi[2] = 360*xi2 - 360*xi + 60;
		phi[3] = -180*xi2 + 168*xi - 24;
		phi[4] = 30*xi2 - 24*xi + 3;
		break;
	case 7: // AF7: uL sL ssL uM uR sR ssR
		phi[0] = 2100*xi4 - 4320*xi3 + 2700*xi2 - 480*xi;
		phi[1] = 420*xi4 - 900*xi3 + 600*xi2 - 120*xi;
		phi[2] = 35*xi4 - 80*xi3 + 60*xi2 - 16*xi + 1;
		phi[3] = -4200*xi4 + 8400*xi3 - 5040*xi2 + 840*xi;
		phi[4] = 2100*xi4 - 4080*xi3 + 2340*xi2 - 360*xi;
		phi[5] = -420*xi4 + 780*xi3 - 420*xi2 + 60*xi;
		phi[6] = 35*xi4 - 60*xi3 + 30*xi2 - 4*xi;
		break;
	}
}
void FindIntBasis(int p, double xi, double *phi)
{
	double xi2 = xi*xi;
	double xi3 = xi*xi*xi;
	double xi4 = xi*xi*xi*xi;
	double xi5 = xi*xi*xi*xi*xi;
	double xi6 = xi*xi*xi*xi*xi*xi;
	double xi7 = xi*xi*xi*xi*xi*xi*xi;
	switch (p)
	{
	case 3: // AF3: uL uM uR
		phi[0] = xi3 - 2*xi2 + xi;
		phi[1] = -2*xi3 + 3*xi2;
		phi[2] = xi3 - xi2;
		break;
	case 5: // AF5: uL sL uM uR sR
		phi[0] = -3*xi5 + 8*xi4 - 6*xi3 + xi;
		phi[1] = -xi5/2 + 3*xi4/2 - 3*xi3/2 + xi2/2;
		phi[2] = 6*xi5 - 15*xi4 + 10*xi3;
		phi[3] = -3*xi5 + 7*xi4 - 4*xi3;
		phi[4] = xi5/2 - xi4 + xi3/2;
		break;
	case 7: // AF7: uL sL ssL uM uR sR ssR
		phi[0] = 10*xi7 - 36*xi6 + 45*xi5 - 20*xi4 + xi;
		phi[1] = 2*xi7 - 15*xi6/2 + 10*xi5 - 5*xi4 + xi2/2;
		phi[2] = xi7/6 - 2*xi6/3 + xi5 - 2*xi4/3 + xi3/6;
		phi[3] = -20*xi7 + 70*xi6 - 84*xi5 + 35*xi4;
		phi[4] = 10*xi7 - 34*xi6 + 39*xi5 - 15*xi4;
		phi[5] = -2*xi7 + 13*xi6/2 - 7*xi5 + 5*xi4/2;
		phi[6] = xi7/6 - xi6/2 + xi5/2 - xi4/6;
		break;
	}
}

double reconstruct(int p, double* phi, CELL cell)
{
	double U = 0.0;
	switch (p)
	{
	case 3:
		U = phi[0] * cell.uL \
		  + phi[1] * cell.uM \
		  + phi[2] * cell.uR;
		break;
	case 5:
		U = phi[0] * cell.uL \
		  + phi[1] * cell.sL \
		  + phi[2] * cell.uM \
		  + phi[3] * cell.uR \
		  + phi[4] * cell.sR;
		break;
	case 7:
		U = phi[0] * cell.uL \
		  + phi[1] * cell.sL \
		  + phi[2] * cell.ssL \
		  + phi[3] * cell.uM \
		  + phi[4] * cell.uR \
		  + phi[5] * cell.sR \
		  + phi[6] * cell.ssR;
		break;
	}
	return U;
}

double FindU(int p, CELL cell, double xi) // corrected
{
	double *phi;
	phi = new double[p];
	FindBasis(p, xi, phi);
	double U = reconstruct(p, phi, cell);
	return U;
}

double FindIntU(int p, CELL cell, double xi)
{
	double *intphi;
	intphi = new double[p];
	FindIntBasis(p, xi, intphi);
	double intU = reconstruct(p, intphi, cell);
	return intU;
}

double FindDiffU(int p, CELL cell, double xi)
{
	double *diffphi;
	diffphi = new double[p];
	FindDiffBasis(p, xi, diffphi);
	double diffU = reconstruct(p, diffphi, cell);
	return diffU;
}

double FindDiff2U(int p, CELL cell, double xi)
{   
	double *diff2phi;
	diff2phi = new double[p];
	FindDiff2Basis(p, xi, diff2phi);
	double diff2U = reconstruct(p, diff2phi, cell);
	return diff2U;
}