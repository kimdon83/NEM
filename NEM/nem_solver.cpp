#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>

#include "nem_solver.h"
//#include "sparce_matrix.h"

#include "define.h"
//Constructor of a nemSolver object
CnemSolver::CnemSolver(void)
{
    //Initialize member variables
    nDim = 1;
    meshSize[XDIR] = 1.;
    meshSize[YDIR] = 1.;
    meshSize[ZDIR] = 1.;
	memset(albedo, 0, sizeof(real)*6);
    ke = 1.e10;
    omega = 1.0;
	method = 0;
    nHigherOrder = 0;

    keff = 1.;
    Flux = NULL;
}

//Destructor of a nemSolver object
CnemSolver::~CnemSolver(void)
{
    //Erase memories
}

//Read an input file
void CnemSolver::ReadInput(char *filename)
{
    char card[CARDLEN];

    ifstream inFp(filename); //Make input file stream

    //Read input
    while(inFp.peek()!=EOF)
    {
        inFp>>card;

        //Read "Options" section
          if(!strcmp(card, "Options")) this->ReadOptions(inFp);
        
        //Read "CXLibrary" section
        else if(!strcmp(card, "CXLibrary")) CXMan.ReadCXTables(inFp);
		//Read "Geometry" section
		else if (!strcmp(card, "Geometry")) Str.ReadGeometry(inFp, nDim);
    }

	//Prepare calculation
	Str.LinkCXTable(&CXMan);
}

//Read the calculation options
void CnemSolver::ReadOptions(istream &ins)
{
    int i;
    char oneChr, card[CARDLEN];
    bool bEnd = false;

    //Read '('
    ins >> oneChr;

    while(bEnd==false)
    {
        ins >> card;

        if(!strcmp(card, "Dimension")) {
            ins >> nDim;
        }
        else if(!strcmp(card,"MeshSize")) {
            for(i=0; i<nDim; i++) ins >> meshSize[i]; //XDIR, YDIR, ZDIR
        }
		else if(!strcmp(card, "Albedo")) {
			for(i=0; i<nDim*2; i++) ins >> albedo[i]; //XL, XR, YL, YR, ZL, ZR
		}
        else if(!strcmp(card, "ke")) {
            ins >> ke;
        }
        else if(!strcmp(card, ENDSTR)) bEnd = true;
    }
}

//Analyze the core problem by the nem method
void CnemSolver::AnalyzeCore(void)
{
	int i;
	int nIt;
	int nMeshAll;
	

	int **nMeshX = NULL, *nMeshXY = NULL;
	//allocate memory for number of mesh
	nMeshXY = new int[Str.nzNode];
	memset(nMeshXY, 0, sizeof(int)*Str.nyNode);


	nMeshX = new int*[Str.nzNode];
	for (i = 0; i < Str.nzNode; i++) {
		nMeshX[i] = new int[Str.nyNode]; 
		memset(nMeshX[i], 0, sizeof(int)*Str.nyNode);
	}
	
	nMeshAll=Str.CalculateNumberOfMeshes(meshSize, nMeshX, nMeshXY);
	Flux.SetDimension(nMeshAll*NUM_GRP);
	//-------------------------------------------------------------------------
	Str.SetCoefficient(meshSize,nMeshAll);
	//-------------------------------------------------------------------------
	// Solve the eigenvalue equation:
	nIt=Str.ExecutePowerIterationMethod(albedo,keff,Flux,nMeshX,nMeshXY);

    //cout<<"Number of Iterations = "<<setw(4)<<nIt<<", Keff = "<<keff<<endl;

    //Erase memories
	if (nMeshX) {
		for (i = 0; i < Str.nzNode; i++)
			delete[] nMeshX[i];
		delete[] nMeshX;
	}
	if (nMeshXY)
		delete[] nMeshXY;
}


