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
		phi[0] = 3 * xi2 - 4 * xi + 1;
		phi[1] = -6 * xi2 + 6 * xi;
		phi[2] = 3 * xi2 - 2 * xi;
		break;
	case 2: // DG2: uL sL uM uR sR
		phi[0] = -15 * xi4 + 32 * xi3 - 18 * xi2 + 1;
		phi[1] = -5 * xi4 / 2 + 6 * xi3 - 9 * xi2 / 2 + xi;
		phi[2] = 30 * xi4 - 60 * xi3 + 30 * xi2;
		phi[3] = -15 * xi4 + 28 * xi3 - 12 * xi2;
		phi[4] = 5 * xi4 / 2 - 4 * xi3 + 3 * xi2 / 2;
		break;
	case 3: // DG3: uL sL ssL uM uR sR ssR
		phi[0] = 70 * xi6 - 216 * xi5 + 225 * xi4 - 80 * xi3 + 1;
		phi[1] = 14 * xi6 - 45 * xi5 + 50 * xi4 - 20 * xi3 + xi;
		phi[2] = 7 * xi6 / 6 - 4 * xi5 + 5 * xi4 - 8 * xi3 / 3 + xi2 / 2;
		phi[3] = -140 * xi6 + 420 * xi5 - 420 * xi4 + 140 * xi3;
		phi[4] = 70 * xi6 - 204 * xi5 + 195 * xi4 - 60 * xi3;
		phi[5] = -14 * xi6 + 39 * xi5 - 35 * xi4 + 10 * xi3;
		phi[6] = 7 * xi6 / 6 - 3 * xi5 + 5 * xi4 / 2 - 2 * xi3 / 3;
		break;
	}
}

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
		psi[1] = 10 * xi2 - 8 * xi + 1;
		break;
	case 2: // DG2: uM uR sR
		psi[0] = 1;
		psi[1] = 126 * xi4 - 224 * xi3 + 126 * xi2 - 24 * xi + 1;
		psi[2] = 105 * xi4 - 196 * xi3 + 231 * xi2 / 2 - 23 * xi + 1;
		break;
	case 3: // DG3: uM uR sR ssR
		psi[0] = 1;
		psi[1] = 1716 * xi6 - 4752 * xi5 + 4950 * xi4 - 2400 * xi3 + 540 * xi2 - 48 * xi + 1;
		psi[2] = 3003 * xi6 / 2 - 4257 * xi5 + 9075 * xi4 / 2 - 2250 * xi3 + 1035 * xi2 / 2 - 47 * xi + 1;
		psi[3] = 4004 * xi6 / 3 - 3850 * xi5 + 4180 * xi4 - 6340 * xi3 / 3 + 496 * xi2 - 46 * xi + 1;
		break;
	}
}

void FindDiffTest(int p, double xi, double *diffpsi)
{
	double xi2 = xi*xi;
	double xi3 = xi*xi*xi;
	double xi4 = xi*xi*xi*xi;
	double xi5 = xi*xi*xi*xi*xi;
	switch (p)
	{
	case 1: // DG1: uM uR
		diffpsi[0] = 0;
		diffpsi[1] = 20 * xi - 8;
		break;
	case 2: // DG2: uM uR sR
		diffpsi[0] = 0;
		diffpsi[1] = 504 * xi3 - 672 * xi2 + 252 * xi - 24;
		diffpsi[2] = 420 * xi3 - 588 * xi2 + 231 * xi - 23;
		break;
	case 3: // DG3: uM uR sR ssR
		diffpsi[0] = 0;
		diffpsi[1] = 10296 * xi5 - 23760 * xi4 + 19800 * xi3 - 7200 * xi2 + 1080 * xi - 48;
		diffpsi[2] = 9009 * xi5 - 21285 * xi4 + 18150 * xi3 - 6750 * xi2 + 1035 * xi - 47;
		diffpsi[3] = 8008 * xi5 - 19250 * xi4 + 16720 * xi3 - 6340 * xi2 + 992 * xi - 46;
		break;
	}
}

double FindU(int p, double *stateLocal, double xi) // corrected
{
	double U = 0.0;
	double *phi;
	phi = new double[2*p+1];
	FindBasis(p, xi, phi);
	for (int i = 0;i < 2*p + 1;i++)
		U += phi[i] * stateLocal[i];
	delete phi;
	return U;
}

