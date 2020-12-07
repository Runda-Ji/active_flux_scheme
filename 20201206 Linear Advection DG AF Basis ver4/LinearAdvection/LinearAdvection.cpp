// LinearAdvection.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "BasisFunc.cpp"
#include "Quadrature.cpp"
#include "BasicMatrixOperation.cpp"
#include "Printing.cpp"

void InitialCondition(MESH *mesh, STATE *state)
{
	int p = mesh->p;
	double Dx = mesh->Dx;
	int NumCells = mesh->NumCells;
	mesh->cells = new CELL[NumCells];

	state->p = p;
	state->NumCells = NumCells;
	state->u = new double[p + (p + 1)*NumCells];

	for (int k = 0; k<NumCells; k++)
	{
		mesh->cells[k].xL = k*Dx;
		mesh->cells[k].xR = (k + 1)*Dx;
		for (int i = 0; i<p + 1; i++)
			state->u[p + k*(p + 1) + i] = 0.0;
	}
	int k = 9;
	state->u[p + k*(p + 1) + 0] = 0.5; // uM
	state->u[p + k*(p + 1) + 1] = 1.0; // uR
	for (k = 10; k < 29; k++)
	{
		state->u[p + k*(p + 1) + 0] = 1.0; // uM
		state->u[p + k*(p + 1) + 1] = 1.0; // uR
	}
	k = 29;
	state->u[p + k*(p + 1) + 0] = 0.5; // uM
									   // periodic boundary condition
	for (int i = 0; i < p; i++)
		state->u[i] = state->u[(NumCells - 1)*(p + 1) + i];
}

double FindInteriorFlux(int p, double *StateLocal, double xi, double a)
{
	double u = FindU(p, StateLocal, xi);
	double f = a*u; // For linear advection
	return f;
}

void InteriorContribution(STATE StateGlobal, double *R, double a)
{
	int p = StateGlobal.p;
	int NumCells = StateGlobal.NumCells;
	QUADRATURE quad;
	SelectQuadPts(p, &quad);
	double *StateLocal;
	StateLocal = new double[2*p + 1];
	// in each iteration, we will first run InteriorContribution and then run EdgeContribution
	// so we need to clean up vector R here, p state for left-most boundary & p+1 states for each cell
	for (int i = 0; i < p + NumCells*(p + 1); i++)
		R[i] = 0;
	
	// loop over all cells
	for (int k=0; k<NumCells; k++) // k is the index of cell
	{
		// state: 2p + 1: uL, sL, ssL, uM, uR, sR, ssR
		for (int j=0; j<2*p+1; j++)
			StateLocal[j] = StateGlobal.u[k*(p+1)+j]; // Inject StateGlobal to StateLocal									  
		// test func: p + 1: psiuM, psiuR, psisR, psissR
		for (int i=0; i<p+1; i++)
		{
			for (int q=0; q<quad.n; q++) // q is the index of quad point
			{
				double* DpsiDxi;
				DpsiDxi = new double[p + 1];
				FindDiffTest(p, quad.x[q], DpsiDxi); // DpsiDxi is a (p+1) vector
				double f = FindInteriorFlux(p, StateLocal, quad.x[q], a); // find the interior flux at quad point q
				R[p+k*(p+1)+i] -= DpsiDxi[i] * f * quad.w[q]; // update the residual for the i^th basis func, sum over q
				delete DpsiDxi;
			}
		}
	}
	delete StateLocal;
}

double FindEdgeFlux(int p, double *StateLocalL, double *StateLocalR, double a)
{   // a > 0, f = a*uL
	// both *StateLocalL and *StateLocalR include 2p + 1 states
	double u = FindU(p, StateLocalL, 1);
	double f = a*u; // For linear advection;
	return f;
}

void EdgeContribution(STATE StateGlobal, double *R, double a)
{
	int p = StateGlobal.p;
	int NumCells = StateGlobal.NumCells;
	// state: 2p + 1: uL, sL, ssL, uM, uR, sR, ssR
	double *StateCurr, *StateNext;
	StateCurr = new double[2*p + 1];
	StateNext = new double[2*p + 1];
	// test func: p + 1: psiuM, psiuR, psisR, psissR
	double *psiCurr, *psiNext;
	psiCurr = new double[p + 1];
	psiNext = new double[p + 1];
	FindTest(p, 1.0, psiCurr);
	FindTest(p, 0.0, psiNext);
	// loop over cells
	for (int k = 0;k < NumCells;k++)
	{
		int kNext = (k + 1) % NumCells;
		// state: 2p + 1: uL, sL, ssL, uM, uR, sR, ssR
		for (int i = 0;i < 2*p + 1;i++)
		{
			StateCurr[i] = StateGlobal.u[    k*(p + 1) + i]; // Inject StateGlobal to current cell
			StateNext[i] = StateGlobal.u[kNext*(p + 1) + i]; // Inject StateGlobal to next cell
		}
		double f = FindEdgeFlux(p, StateCurr, StateNext, a);
		// residual: p + 1 only update uM, uR, sR, ssR
		for (int i = 0;i < p + 1;i++)
		{
			R[p +     k*(p + 1) + i] += psiCurr[i] * f;
			R[p + kNext*(p + 1) + i] -= psiNext[i] * f;
		}
	}
	delete StateCurr;
	delete StateNext;
	delete psiCurr;
	delete psiNext;
}

void CopyState(STATE stateIn, STATE *stateOut) {
	int p = stateIn.p;
	int NumCells = stateIn.NumCells;
	for (int i = 0; i < p + NumCells*(p + 1); i++)
		stateOut->u[i] = stateIn.u[i];
}


void SubStepRK4(STATE State0, double a, double coef, double *iM, double *ResGlobal, STATE *StateOut)
{
	int p = State0.p;
	int NumCells = State0.NumCells;
	StateOut->p = p;
	StateOut->NumCells = NumCells;
	// only need to update uM, uR, sR, ssR
	double *ResLocal, *StateLocal;
	ResLocal   = new double[p + 1];
	StateLocal = new double[p + 1];
	// loop over cells
	for (int k = 0;k < NumCells;k++)
	{
		for (int i = 0; i < p + 1; i++)
		{
			ResLocal[i]   = ResGlobal[p + k *(p + 1) + i];
			StateLocal[i] = State0.u[p + k *(p + 1) + i];
		}
		// StateLocal = StateLocal + coef * iM * ResLocal;
		MVProd(coef, iM, ResLocal, StateLocal, p + 1, p + 1);
		for (int i = 0; i < p + 1; i++)
			StateOut->u[p + k *(p + 1) + i] = StateLocal[i];
	}
	// left-most boundary
	for (int i = 0; i < p; i++)
		StateOut->u[i] = StateOut->u[p + (NumCells - 1)*(p + 1) + i];
	delete ResLocal;
	delete StateLocal;
}


void DiscontinuousGalerkin(MESH *mesh, STATE *stateExternal, double a, double t, double CFL)
{	
	int p = mesh->p;
	int NumCells = mesh->NumCells;
	double Dx = mesh->Dx;
	double Dt = CFL*Dx / a; // time step
	
	int NTimeStep = int(t / Dt); // # of time step
	// int NTimeStep = 100000; // debug

	double *M, *iM;
	M = new double[(p + 1)*(p + 1)];
	iM = new double[(p + 1)*(p + 1)];
	MassMatrix(p, Dx, M);
	// PrintMatrix(M, p + 1, p + 1);
	Inverse(M, iM, p + 1);

	double *R;
	R = new double[p + NumCells*(p + 1)];

	STATE state[5];
	for (int i = 0; i < 5; i++)
	{
		state[i].p = p;
		state[i].NumCells = NumCells;
		state[i].u = new double[p + NumCells*(p + 1)];
	}
		
	CopyState(*stateExternal, &state[0]);
	
	for (int i = 0;i < NTimeStep;i++)
	{
		// PrintState(state[0]);
		// step 1
		InteriorContribution(state[0], R, a);		
		EdgeContribution(state[0], R, a);
		// double ResNorm = ComputeNorm("1", R, 1, NumCells*(p + 1)); // compute residual
		// printf("iter=%d res=%.16e\n", i, ResNorm);                 // print residual
		SubStepRK4(state[0], a, -Dt / 4.0, iM, R, &state[1]);      // compute state[1]
		// step 2
		InteriorContribution(state[1], R, a);
		EdgeContribution(state[1], R, a);
		SubStepRK4(state[0], a, -Dt / 3.0, iM, R, &state[2]);      // compute state[2]
		// step 3
		InteriorContribution(state[2], R, a);
		EdgeContribution(state[2], R, a);
		SubStepRK4(state[0], a, -Dt / 2.0, iM, R, &state[3]);      // compute state[3]
		// step 4
		InteriorContribution(state[3], R, a);
		EdgeContribution(state[3], R, a);
		SubStepRK4(state[0], a, -Dt / 1.0, iM, R, &state[4]);      // compute state[4]
		CopyState(state[4], &state[0]);                            // update state[0]
	}
	CopyState(state[0], stateExternal);
	delete R;
	for (int i = 0;i < 5;i++)
		delete state[i].u;
}

int main()
{
	MESH mesh;
	STATE state;
	mesh.NumCells = 40;
	mesh.Dx = 1;

	int pList[] = { 1, 2, 3 };
	double CFLMaxList[] = {0.464, 0.235, 0.145};
	double coefList[] = {0.25, 0.5, 0.75};

	for (int i = 0; i < 3; i++)
	{
		mesh.p = pList[i];
		InitialCondition(&mesh, &state);
		double a = 1.0;
		double t = mesh.NumCells*mesh.Dx / a;
		double CFLMax = CFLMaxList[i];	
		for (int j = 0; j < 3; j++)
		{
			double CFL = coefList[j] * CFLMax;
			DiscontinuousGalerkin(&mesh, &state, a, t, CFL);
			// PrintState(state);
			SolutionOutput(state, "DG" + to_string(state.p) + "_CFL" + to_string(CFL));
		}		
	}
	system("pause");
	return 0;
}

