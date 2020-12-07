// LinearAdvection.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "BasisFunc.cpp"

double CharacterOri(double a, double Dt, double Dx)
{
	double xi0 = 1 - a * (Dt / Dx); // 1 - nu
	return xi0;
} 

double AvgFlux(int p, CELL cell, double a, double Dt, double Dx)
{
	double xi1 = 1;
	double xi0 = CharacterOri(a,Dt,Dx);
	double flux = (Dx / Dt) * (FindIntU(p, cell, xi1) - FindIntU(p, cell, xi0));
	return flux;
}

double UpdateCellAvg(CELL cell, double Dt, double Dx, double fR, double fL)
{
	double uM = cell.uM - (Dt / Dx)*(fR - fL);
	return uM;
}

void UpdateRightEdge(int p, STATE *NewState, CELL cell, double a, double Dt, double Dx)
{
	double xi1 = 1;
	double xi0 = CharacterOri(a, Dt, Dx);
	NewState->u = FindU(p, cell, xi0);
	NewState->s = FindDiffU(p, cell, xi0);
	NewState->ss = FindDiff2U(p, cell, xi0);
}

void InitialConditionStep(MESH *mesh)
{
	// set up coordinates, domain [0,40]
	mesh->Dx = 1;
	for (int i = 0;i < mesh->NumCells;i++)
	{
		mesh->cells[i].xL = 0 + i * mesh->Dx;
		mesh->cells[i].xR = 0 + (i + 1) * mesh->Dx;
	}
	// projecting ini conditions
	int tri1 = 9;
	int tri2 = 29;
	for (int i = 0; i < tri1; i++)
	{
		mesh->cells[i].uL = 0;
		mesh->cells[i].uR = 0;
		mesh->cells[i].uM = 0;
	}
	mesh->cells[tri1].uL = 0;
	mesh->cells[tri1].uR = 1;
	mesh->cells[tri1].uM = 0.5;
	for (int i = tri1+1; i < tri2; i++)
	{
		mesh->cells[i].uL = 1;
		mesh->cells[i].uR = 1;
		mesh->cells[i].uM = 1;
	}
	mesh->cells[tri2].uL = 1;
	mesh->cells[tri2].uR = 0;
	mesh->cells[tri2].uM = 0.5;
	for (int i = tri2 + 1; i < mesh->NumCells; i++)
	{
		mesh->cells[i].uL = 0;
		mesh->cells[i].uR = 0;
		mesh->cells[i].uM = 0;
	}
	for (int i = 0;i < mesh->NumCells;i++)
	{
		mesh->cells[i].sL = 0;
		mesh->cells[i].sR = 0;
		mesh->cells[i].ssL = 0;
		mesh->cells[i].ssR = 0;
	}
}

void SolutionOutput(MESH mesh, string title)
{
	ofstream output;
	string fname = "data/" + title + ".txt";
	output.open(fname);
	for (int i = 0;i < mesh.NumCells;i++)
	{
		/*
		double xi = 0.5;
		double u = FindU(mesh.p,mesh.cells[i], xi);
		output << u << endl;
		*/
		double uM = mesh.cells[i].uM;
		output << uM << endl;
		
	}
	output.close();
}

void LinearAdvection(MESH *mesh, double a, double t, double CFL)
{
	int p = mesh->p;
	int NumCells = mesh->NumCells;
	double Dx = mesh->Dx;
	double Dt = CFL*Dx / a;
	int NTimeStep = int(t / Dt);
	for (int i = 0;i < NTimeStep;i++) // loop for time marching
	{
		double *NewuM;
		NewuM = new double[NumCells];
		STATE *NewState;
		NewState = new STATE[NumCells];

		for (int j = 0;j < NumCells;j++) // loop over all cells, and find the new value
		{
			CELL cell     = mesh->cells[j];
			CELL cellprev = mesh->cells[(j+NumCells-1) % NumCells];
			double fR = AvgFlux(p, cell, a, Dt, Dx);
			double fL = AvgFlux(p, cellprev, a, Dt, Dx);
			NewuM[j] = UpdateCellAvg(cell, Dt, Dx, fR, fL);
			UpdateRightEdge(p, &(NewState[j]), cell, a, Dt, Dx);
		}
		for (int j = 0;j < NumCells;j++) // loop over all cells, and actually update the cells
		{
			int jprev = (j+NumCells-1) % NumCells;
			mesh->cells[j].uM  = NewuM[j];
			mesh->cells[j].uR  = NewState[j].u;
			mesh->cells[j].sR  = NewState[j].s;
			mesh->cells[j].ssR = NewState[j].ss;
			mesh->cells[j].uL  = NewState[jprev].u;
			mesh->cells[j].sL  = NewState[jprev].s;
			mesh->cells[j].ssL = NewState[jprev].ss;
		}
	}
}

int main()
{
	CELL cell;
	MESH mesh;
	mesh.NumCells = 40;
	mesh.cells = new CELL[mesh.NumCells];

	int pList[] = { 3, 5, 7 };
	double CFLList[] = { 0.25, 0.5, 0.75 };

	for (int i = 0; i < 3; i++)
	{
		mesh.p = pList[i];
		InitialConditionStep(&mesh);
		double a = 1.0; // wave speed
		double t = mesh.NumCells*mesh.Dx / a;
		for (int j = 0; j < 3; j++)
		{
			double CFL = CFLList[j];
			LinearAdvection(&mesh, a, t, CFL);
			SolutionOutput(mesh, "AF" + to_string(mesh.p) + "_CFL" + to_string(CFL));
		}
	}
	
	
	system("pause");
    return 0;
}

