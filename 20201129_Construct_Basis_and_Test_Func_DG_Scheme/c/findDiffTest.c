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
		diffpsi[1] = 20*xi - 8;
		break;
	case 2: // DG2: uM uR sR
		diffpsi[0] = 0;
		diffpsi[1] = 504*xi3 - 672*xi2 + 252*xi - 24;
		diffpsi[2] = 420*xi3 - 588*xi2 + 231*xi - 23;
		break;
	case 3: // DG3: uM uR sR ssR
		diffpsi[0] = 0;
		diffpsi[1] = 10296*xi5 - 23760*xi4 + 19800*xi3 - 7200*xi2 + 1080*xi - 48;
		diffpsi[2] = 9009*xi5 - 21285*xi4 + 18150*xi3 - 6750*xi2 + 1035*xi - 47;
		diffpsi[3] = 8008*xi5 - 19250*xi4 + 16720*xi3 - 6340*xi2 + 992*xi - 46;
		break;
	}
}