void FindBasis(int p, double xi, double *phi)
{
	double xi2 = xi*xi;
	double xi3 = xi*xi*xi;
	double xi4 = xi*xi*xi*xi;
	double xi5 = xi*xi*xi*xi*xi;
	double xi6 = xi*xi*xi*xi*xi*xi;
	switch (p)
	{
	case 1: // DG1: uL uM uR
		phi[0] = 3*xi2 - 4*xi + 1;
		phi[1] = -6*xi2 + 6*xi;
		phi[2] = 3*xi2 - 2*xi;
		break;
	case 2: // DG2: uL sL uM uR sR
		phi[0] = -15*xi4 + 32*xi3 - 18*xi2 + 1;
		phi[1] = -5*xi4/2 + 6*xi3 - 9*xi2/2 + xi;
		phi[2] = 30*xi4 - 60*xi3 + 30*xi2;
		phi[3] = -15*xi4 + 28*xi3 - 12*xi2;
		phi[4] = 5*xi4/2 - 4*xi3 + 3*xi2/2;
		break;
	case 3: // DG3: uL sL ssL uM uR sR ssR
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