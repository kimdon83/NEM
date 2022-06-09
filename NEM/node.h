#ifndef _CELL_H_
#define _CELL_H_

#include "define.h"
#include "cxman.h"
#include "sparce_matrix.h"

/**************************************************************************************************
 * Unit Node
**************************************************************************************************/
class CNode {
private:
    CNode *NeighborNode[NEIGH_NUM]; //[NEIGHBOR_INDEX]
	int nMesh[3]; //XDIR, YDIR, ZDIR

	real meshVol;
    real meshArea[3];
public:
	real meshSize[3];

	real OutgoingCurrent[6][2];
	real InComingCurrent[6][2];
	real L[3][2];	// transverse leakeage direction ans energy group
	real Qvalue[4][3][2]; // index and direction and energy group
	real Cvalue[5][3][2]; // index and direction and energy group
	real Kvalue[3][3][2]; // index and direction and energy group
	real Flux[2];
	real Nmat[3]; // N11 N21 N22 .   N12는 원소가 없으므로 
	real Nmat2[2][2];
	real Mmat1[3][2];
	real Mmat2[3][2];
	real beta[3][2]; //beta_ug ; 
	real Qmat[2]; // energy group
	real Smat[2]; // energy group

	real source[2]; // energy group

public:
    int cxtblId; //xs table id given in the input file
	int configId; // 
    CCXTable *CX; // This Cell's cross-sections
	//CCXTable **CXs; // adresses of the cross sections of Node

public:
    CNode(void); //Constructor
	~CNode(void) { delete[] CX; delete[] NeighborNode; } //Destructor

	void SetMeshSize(real size[3]) {
		memcpy(meshSize, size, sizeof(real)*3); 
		meshVol = meshSize[XDIR]*meshSize[YDIR]*meshSize[ZDIR];
        meshArea[XDIR] = meshSize[YDIR]*meshSize[ZDIR];
        meshArea[YDIR] = meshSize[XDIR]*meshSize[ZDIR];
        meshArea[ZDIR] = meshSize[XDIR]*meshSize[YDIR];
	}
	real GetMeshSize(int dir) { return meshSize[dir]; }

	real GetMeshVol(void) { return meshVol; }

	void SetMeshNumbers(int nX, int nY, int nZ) { 
		nMesh[XDIR] = nX; nMesh[YDIR] = nY; nMesh[ZDIR] = nZ;
	}
	int GetNumberOfMesh(int dir) { return nMesh[dir]; }
	int GetNumberOfMeshes(void) { return nMesh[XDIR]*nMesh[YDIR]*nMesh[ZDIR]; }

	int GetCXTableID(void) { return cxtblId; }
	void SetCXTable(CCXTable *pCX) { CX = pCX; }
	void SetCXTables(CCXTable *pCX) {  }

	void LinkNodeXL(CNode *pNode) { NeighborNode[XL] = pNode; }
	void LinkNodeXR(CNode *pNode) { NeighborNode[XR] = pNode; }
	void LinkNodeYL(CNode *pNode) { NeighborNode[YL] = pNode; }
	void LinkNodeYR(CNode *pNode) { NeighborNode[YR] = pNode; }
	void LinkNodeZL(CNode *pNode) { NeighborNode[ZL] = pNode; }
	void LinkNodeZR(CNode *pNode) { NeighborNode[ZR] = pNode; }

	void SetCoefficient(void);
	void Nem(real *Albedo, real keff);

};

/**************************************************************************************************
 * Structure 
**************************************************************************************************/
class CStructure {
public:
    CNode ****Node; //[Z][Y][X]
    int nyNode;
    int nxNode;
	int nzNode;

    real *xNodeSize;
    real *yNodeSize;
	real *zNodeSize;

	int **nodeTypeId;
	int nNodeType;

public:
    CStructure(void); //Constructor
    ~CStructure(void); //Destructor
	

	//Link CXTable into each Node
	void LinkNeighborNodes(void);
	void LinkCXTable(CCXManager *pCXMan);

    CNode* GetNode(int zId,int yId, int xId) { return Node[zId][yId][xId]; }

public:
    void ReadGeometry(istream &ins, int nDim); //Read "Geometry" section in a input file

	//Return the dimension of the system matrix
	int CalculateNumberOfMeshes(real *MeshSize, int **nMeshX, int *nMeshXY);
	void SetCoefficient(real * MeshSize,int nMeshAll);
	int ExecutePowerIterationMethod(real *albedo,real keff,CVector &FluxOut
		, int ** nMeshX, int * nMeshXY);
	

};
#endif