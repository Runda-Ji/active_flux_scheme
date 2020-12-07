void FindTest(int p, double xi, double *psi)
{
	double xi2 = xi*xi;
	double xi3 = xi*xi*xi;
	double xi4 = xi*xi*xi*xi;
	double xi5 = xi*xi*xi*xi*xi;
	double xi6 = xi*xi*xi*xi*xi*xi;
	switch (p)
	{
	case 1: // DG1: uM uR
		psi[0] = 1;
		psi[1] = 10*xi2 - 8*xi + 1;
		break;
	case 2: // DG2: uM uR sR
		psi[0] = 1;
		psi[1] = 126*xi4 - 224*xi3 + 126*xi2 - 24*xi + 1;
		psi[2] = 105*xi4 - 196*xi3 + 231*xi2/2 - 23*xi + 1;
		break;
	case 3: // DG3: uM uR sR ssR
		psi[0] = 1;
		psi[1] = 1716*xi6 - 4752*xi5 + 4950*xi4 - 2400*xi3 + 540*xi2 - 48*xi + 1;
		psi[2] = 3003*xi6/2 - 4257*xi5 + 9075*xi4/2 - 2250*xi3 + 1035*xi2/2 - 47*xi + 1;
		psi[3] = 4004*xi6/3 - 3850*xi5 + 4180*xi4 - 6340*xi3/3 + 496*xi2 - 46*xi + 1;
		break;
	}
}