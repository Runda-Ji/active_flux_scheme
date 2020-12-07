void Inverse(double*A, double*A_inv, int N)
{
	int L = N*N;
	int *IPIV = new int[N + 1];
	double *WORK = new double[L];
	int INFO;

	dlacpy_("", &N, &N, A, &N, A_inv, &N);
	dgetrf_(&N, &N, A_inv, &N, IPIV, &INFO);
	dgetri_(&N, A_inv, &N, IPIV, WORK, &L, &INFO);

	delete IPIV;
	delete WORK;
}

void MMProd(double*A, double*B, double*C, int m, int k, int n)
{
	// A = m*k, B = k*n, C = m*n
	// C = alpha*A*B + beta*C 
	double alpha = 1;
	double beta = 0;
	dgemm_("N", "N", &m, &n, &k, &alpha, A, &m, B, &k, &beta, C, &n);
}

void MVProd(double alpha, double*A, double*x, double*y, int m, int n)
{
	// A = m*n, x = m*1, y = m*1
	// y = alpha*A*x + beta*y
	double beta = 1;
	int incx = 1;
	int incy = 1;
	dgemv_("N", &m, &n, &alpha, A, &m, x, &incx, &beta, y, &incy);
}

double ComputeNorm(char *c, double *A, int m, int n)
{
	double norm;
	double *work;
	work = new double[m];
	norm = dlange_(c, &m, &n, A, &m, work);
	delete work;
	return norm;
}

void MassMatrix(int p, double Dx, double*M)
{
	QUADRATURE quad;
	SelectQuadPts(p, &quad);
	double *psi, *phi;
	psi = new double[p + 1];   // test func: uM uR sR ssR
	phi = new double[2*p + 1]; // basis func: uL sL ssL uM uR sR ssR
	for (int i=0; i<p+1; i++)
	{
		for (int j=0; j<p+1; j++)
		{
			M[i*(p + 1) + j] = 0;
			for (int q=0; q<quad.n; q++)
			{
				FindTest(p, quad.x[q], psi);
				FindBasis(p, quad.x[q], phi);
				M[i*(p + 1) + j] += psi[i] * phi[p+j] * quad.w[q];
			}
			M[i*(p + 1) + j] *= Dx;
		}
	}
	delete psi;
	delete phi;
}