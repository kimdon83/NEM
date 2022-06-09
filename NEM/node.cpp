#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>

#include "node.h"
#include "vector.h"
#include "sparce_matrix.h"

#include "define.h"
//Constructor of Cell
CNode::CNode(void)
{
	int i;

	//Initialize member variables
	for (i = 0; i < NEIGH_NUM; i++) NeighborNode[i] = NULL;
	for (i = 0; i < 3; i++) {
		nMesh[i] = 1; meshSize[i] = 1.;
		meshArea[i] = 1.;
	}
	for (i = 0; i < 6; i++) {
		OutgoingCurrent[i][0] = 1.;
		InComingCurrent[i][0] = 1.;
		OutgoingCurrent[i][1] = 1.;
		InComingCurrent[i][1] = 1.;
	}	//Outgoing and incoming Current
	for (i = 0; i < 3; i++) {
		L[i][0] = 1.;
		L[i][1] = 1.;
	}	//transverse leakeage

	Flux[0] = 1.;
	Flux[1] = 1.;

	Qmat[0] = 0;
	Qmat[1] = 0;

	Smat[0] = 0;
	Smat[1] = 0;

	meshVol = 1.;
	cxtblId = -1;
	CX = NULL;
}
void CNode::SetCoefficient(void)
{
	int g, i;
	for (g = 0; g < NUM_GRP; g++) {
		for (i = 0; i < 3; i++) {
			beta[i][g] = CX->GetDiffCoeff(g) / meshSize[i];
			Qvalue[0][i][g] = beta[i][g] / (1 + 12 * beta[i][g]);
			Qvalue[1][i][g] = beta[i][g] / (1 + 4 * beta[i][g]) / (1 + 12 * beta[i][g]);
			Qvalue[2][i][g] = (1 - 48 * beta[i][g] * beta[i][g]) / (1 + 4 * beta[i][g]) / (1 + 12 * beta[i][g]);
			Qvalue[3][i][g] = beta[i][g] / (1 + 4 * beta[i][g]);
		}
	}
	for (i = 0; i < 3; i++) {
		Qmat[0] += Qvalue[0][i][0] / meshSize[i];
		Qmat[1] += Qvalue[0][i][1] / meshSize[i];
	}

	//Nmat[3]을 미리 계산한다.
	Nmat[0] = CX->GetCX_ABS(0) + CX->GetScaDiff(0, 1);
	Nmat[1] = -CX->GetScaDiff(0, 1);

	Nmat[2] = CX->GetCX_ABS(1);


}
void CNode::Nem(real *Albedo, real keff)
{
	string abcc;
	int i, g;
	CNode *pNode;
	real alb;
	real surL[6][2];	//L_ulg, L_urg
	real denom;


	//주변의 셀들의 outgoing current로 부터 incoming current를 받아온다.
	for (g = 0; g < NUM_GRP; g++) {
		for (i = 0; i < 6; i++) {
			pNode = NeighborNode[i];
			alb = Albedo[i];
			if (pNode) {
				if ((i == XL) || (i == YL) || (i == ZL)) {
					InComingCurrent[i][g] = pNode->OutgoingCurrent[i + 1][g];
				}
				else {
					InComingCurrent[i][g] = pNode->OutgoingCurrent[i - 1][g];
				}
			}
			else {		//예외처리를 통해 경계조건을 적용한다.
				InComingCurrent[i][g] = OutgoingCurrent[i][g] * (1. - 2. * alb) / (1. + 2. * alb);
				//if (alb == 0) {
				//	InComingCurrent[i][g] = OutgoingCurrent[i][g];
				//}
				//else {
				//	InComingCurrent[i][g] = OutgoingCurrent[i][g] * (1 - 2 * alb) / (1 + 2 * alb);
				//}
			}
		}
	}

	//average transverse leakeage를 계산하여서 L_ulg,L_urg를 계산한다.
	for (g = 0; g < NUM_GRP; g++) {
		L[XDIR][g] = (1. / CX->GetDiffCoeff(g))*
			((OutgoingCurrent[YR][g] + OutgoingCurrent[YL][g] - InComingCurrent[YR][g] - InComingCurrent[YL][g]) / meshSize[YDIR]
				+ (OutgoingCurrent[ZR][g] + OutgoingCurrent[ZL][g] - InComingCurrent[ZR][g] - InComingCurrent[ZL][g]) / meshSize[ZDIR]);
		L[YDIR][g] = (1. / CX->GetDiffCoeff(g))*
			((OutgoingCurrent[XR][g] + OutgoingCurrent[XL][g] - InComingCurrent[XR][g] - InComingCurrent[XL][g]) / meshSize[XDIR]
				+ (OutgoingCurrent[ZR][g] + OutgoingCurrent[ZL][g] - InComingCurrent[ZR][g] - InComingCurrent[ZL][g]) / meshSize[ZDIR]);
		L[ZDIR][g] = (1. / CX->GetDiffCoeff(g))*
			((OutgoingCurrent[XR][g] + OutgoingCurrent[XL][g] - InComingCurrent[XR][g] - InComingCurrent[XL][g]) / meshSize[XDIR]
				+ (OutgoingCurrent[YR][g] + OutgoingCurrent[YL][g] - InComingCurrent[YR][g] - InComingCurrent[YL][g]) / meshSize[YDIR]);
	}
	for (g = 0; g < NUM_GRP; g++) {
		for (i = 0; i < 6; i++) {
			pNode = NeighborNode[i];
			if (pNode) {
				surL[i][g] = (beta[int(i / 2)][g]*L[int(i / 2)][g]+ (pNode-> beta[int(i / 2)][g]*pNode->L[int(i / 2)][g]))/ (pNode->beta[int(i / 2)][g] + beta[int(i / 2)][g]);
			}
			else {
				if (Albedo[i] == 0) {	//reflective B.C
					surL[i][g] = L[int(i / 2)][g];
				}
				else { 			//zero incoming
					//surL[i][g] =0;
					surL[i][g] = L[int(i / 2)][g]/2.;

				}
			}
		}
	}

	//이를 통해 Kvalue를 계산한다.
	for (g = 0; g < NUM_GRP; g++) {
		for (i = 0; i < 3; i++) {
			Kvalue[0][i][g] = L[i][g];
			Kvalue[1][i][g] = (-surL[i * 2][g] + surL[i * 2 + 1][g]) / 2.;
			Kvalue[2][i][g] = (surL[i * 2][g] + surL[i * 2 + 1][g]) / 2. - L[i][g];
		}
	}

	//Cvalue를 계산한다.
	for (g = 0; g < NUM_GRP; g++) {
		for (i = 0; i < 3; i++) {
			Cvalue[0][i][g] = Flux[g];
			Cvalue[1][i][g] = (OutgoingCurrent[i * 2 + 1][g] + InComingCurrent[i * 2 + 1][g])- (OutgoingCurrent[i * 2][g] + InComingCurrent[i * 2][g]);
			Cvalue[2][i][g] = (OutgoingCurrent[i * 2 + 1][g] + InComingCurrent[i * 2 + 1][g])+ (OutgoingCurrent[i * 2][g] + InComingCurrent[i * 2][g]) - Flux[g];
		}
	}
	//계산해놓은 Nmat를 사용하여 이 노드에서 사용할 Nmat2[][]를 계산한다.
	Nmat2[0][0] = Nmat[0] - CX->GetCX_NUFIS(0) / keff;
	Nmat2[0][1] = -CX->GetCX_NUFIS(1) / keff;
	Nmat2[1][0] = Nmat[1];
	Nmat2[1][1] = Nmat[2];
	//Mmat1[3][2], Mmat2[3][2]를 계산한다.
	for (i = 0; i < 3; i++) {
		Mmat1[i][0] = CX->GetDiffCoeff(0)*Kvalue[1][i][0] + Nmat2[0][0] * Cvalue[1][i][0] + Nmat2[0][1] * Cvalue[1][i][1];
		Mmat1[i][1] = CX->GetDiffCoeff(1)*Kvalue[1][i][1] + Nmat2[1][0] * Cvalue[1][i][0] + Nmat2[1][1] * Cvalue[1][i][1];
		Mmat2[i][0] = CX->GetDiffCoeff(0)*Kvalue[2][i][0] + Nmat2[0][0] * Cvalue[2][i][0] + Nmat2[0][1] * Cvalue[2][i][1];
		Mmat2[i][1] = CX->GetDiffCoeff(1)*Kvalue[2][i][1] + Nmat2[1][0] * Cvalue[2][i][0] + Nmat2[1][1] * Cvalue[2][i][1];
	}
	//Cvalue 3,4를 계산한다.
	for (i = 0; i < 3; i++) {
		denom = ((Nmat2[0][0] + 60. * beta[i][0] / meshSize[i])*(Nmat2[1][1] + 60. * beta[i][1] / meshSize[i])-Nmat2[1][0]* Nmat2[0][1]);
		Cvalue[3][i][0] = 10.*((Nmat2[1][1] + 60. * beta[i][1] / meshSize[i])*Mmat1[i][0] - Nmat2[0][1] * Mmat1[i][1]) / denom;
		Cvalue[3][i][1] = 10.*(-Nmat2[1][0] * Mmat1[i][0] + (Nmat2[0][0] + 60. * beta[i][0] / meshSize[i]) * Mmat1[i][1]) / denom;
		denom = ((Nmat2[0][0] + 140. * beta[i][0] / meshSize[i])*(Nmat2[1][1] + 140. * beta[i][1] / meshSize[i])- Nmat2[1][0] * Nmat2[0][1] );
		Cvalue[4][i][0] = 14.*((Nmat2[1][1] + 140. * beta[i][1] / meshSize[i])*Mmat2[i][0] - Nmat2[0][1] * Mmat2[i][1]) / denom;
		Cvalue[4][i][1] = 14.*(-Nmat2[1][0] * Mmat2[i][0] + (Nmat2[0][0] + 140. * beta[i][0] / meshSize[i]) * Mmat2[i][1]) / denom;
	}



	//incoming Current를 이용하여 phi=inv(N+12Q)S를 통해 flux를 계산한다. 
		//여기서 Q는 Qvalue가 아닌 Qmat
		//여기서 N은 매번 바뀌는 Nmat2
		//Smat를 먼저 계산한다.
	Smat[0] = 0.;
	Smat[1] = 0.;
	for (g = 0; g < NUM_GRP;g ++) {
		for (i = 0; i < 3; i++) {
			Smat[g] += (2*Qvalue[0][i][g]*Cvalue[4][i][g]+(InComingCurrent[i*2+1][g]+InComingCurrent[i*2][g])*(1+8*Qvalue[1][i][g]-Qvalue[2][i][g]) ) / meshSize[i];
		}
	}

	denom = (Nmat2[0][0]+12.*Qmat[0])*(Nmat2[1][1] + 12. * Qmat[1])-Nmat2[0][1]*Nmat2[1][0];

	Flux[0] = ( (Nmat2[1][1]+12.*Qmat[1])*Smat[0]-Nmat2[0][1]*Smat[1]) /denom;
	Flux[1] = (-Nmat2[1][0] * Smat[0]+(Nmat2[0][0] + 12. * Qmat[0])*Smat[1]) / denom;

	//flux로 부터 outgoing current를 계산한다. (eq. 6.b-1, 6.b-2)
	for (g = 0; g < NUM_GRP; g++) {
		for (i = 0; i < 3; i++) {
			OutgoingCurrent[i * 2 + 1][g]	= (6. * Flux[g] - Cvalue[4][i][g])*Qvalue[0][i][g] - 8. * InComingCurrent[i * 2][g] * Qvalue[1][i][g]+ InComingCurrent[i * 2 + 1][g] * Qvalue[2][i][g] - Cvalue[3][i][g] * Qvalue[3][i][g];
			OutgoingCurrent[i * 2][g]		= (6. * Flux[g] - Cvalue[4][i][g])*Qvalue[0][i][g] - 8. * InComingCurrent[i * 2 + 1][g] * Qvalue[1][i][g]+ InComingCurrent[i * 2][g] * Qvalue[2][i][g] + Cvalue[3][i][g] * Qvalue[3][i][g];
		}
	}
	
	//////Cvalue를 다시 계산한다.

	source[0] = CX->GetCX_CHI(0)*(CX->GetCX_NUFIS(0)*Flux[0] + CX->GetCX_NUFIS(1)*Flux[1]);
	source[1] = CX->GetCX_CHI(1)*(CX->GetCX_NUFIS(0)*Flux[0] + CX->GetCX_NUFIS(1)*Flux[1]);

	//source[0] = Smat[0];
	//source[1] = Smat[1];
}


//Constructor of the CStructure class
CStructure::CStructure(void)
{
	//Initialize
	Node = NULL;
	nyNode = 0;
	nxNode = 0;
	xNodeSize = NULL;
	yNodeSize = NULL;
}

//Destructor of the CStructure class
CStructure::~CStructure(void)
{
	int i, j, k;

	//Erase memories
	if (Node) {
		for (k = 0; k < nzNode; k++) {
			for (i = 0; i < nyNode; i++) {
				for (j = 0; j < nxNode; j++) delete Node[k][i][j];
				delete[] Node[k][i];
			}
			delete[] Node[k];
		}
		delete[] Node;
	}
	if (nodeTypeId) {
		for (k = 0; k < nNodeType; k++) {
			delete[] nodeTypeId[k];
		}
		delete[] nodeTypeId;
	}
	if (xNodeSize) delete[] xNodeSize;
	if (yNodeSize) delete[] yNodeSize;
	if (zNodeSize) delete[] zNodeSize;

}

//Read "Geometry" section in a input file
void CStructure::ReadGeometry(istream &ins, int nDim)
{
	int i, j, k;
	char oneChr, card[CARDLEN];
	bool bEnd = false;

	//Read '('
	ins >> oneChr;

	while (bEnd == false)
	{
		ins >> card;

		if (!strcmp(card, "NodeNum"))
		{
			if (nDim == 1) {
				nyNode = 1;
				yNodeSize = new real[nyNode];
				yNodeSize[0] = 1.;
			}
			else if (nDim == 2) {

				ins >> nxNode;
				ins >> nyNode;
			}
			else {

				ins >> nxNode;
				ins >> nyNode;
				ins >> nzNode;
			}
			Node = new CNode ***[nzNode];
			for (j = 0; j < nzNode; j++) {
				Node[j] = new CNode**[nyNode];
				for (i = 0; i < nyNode; i++) {
					Node[j][i] = new CNode*[nxNode];
				}
			}

		}
		else if (!strcmp(card, "xNodeSize")) //Read "Size" card
		{
			xNodeSize = new real[nxNode];
			for (i = 0; i < nxNode; i++) ins >> xNodeSize[i];
		}
		else if (!strcmp(card, "yNodeSize")) //Read "Size" card
		{
			yNodeSize = new real[nyNode];
			for (i = 0; i < nyNode; i++) ins >> yNodeSize[i];
		}
		else if (!strcmp(card, "zNodeSize")) //Read "Size" card
		{
			zNodeSize = new real[nzNode];
			for (i = 0; i < nzNode; i++) ins >> zNodeSize[i];
		}



		else if (!strcmp(card, "NodeType"))
		{
			ins >> nNodeType;
			nodeTypeId = new int*[nNodeType];
			for (i = 0; i < nNodeType; i++) {
				nodeTypeId[i] = new int[nzNode];
			}

			for (i = 0; i < nNodeType; i++)
				for (j = 0; j < nzNode; j++)
				{
					ins >> card;
					nodeTypeId[i][j] = atoi(card);
				}
		}
		else if (!strcmp(card, "Configuration"))
		{
			//
			for (i = 0; i < nyNode; i++)
			{
				for (j = 0; j < nxNode; j++)
				{
					ins >> card;
					for (k = 0; k < nzNode; k++) {
						if (strcmp(card, STR_EMPTY)) {

							Node[k][i][j] = new CNode;
							Node[k][i][j]->cxtblId = nodeTypeId[atoi(card)][k];
							//Node[i][j]->configId = atoi(card);
							//Node[i][j]->cxtblId = atoi(card);
						}
						else {
							Node[k][i][j] = NULL;
						}
					}
				}
			}
			//}
		}
		else if (!strcmp(card, ENDSTR)) bEnd = true; //Read ");"
	}

	//Link neighbor nodes
	LinkNeighborNodes();
}

void CStructure::LinkNeighborNodes(void)
{
	int i, j, k;

	for (k = 0; k < nzNode; k++) {
		for (i = 0; i < nyNode; i++) {
			for (j = 0; j < nxNode; j++) {
				if (Node[k][i][j]) {
					//X-LEFT
					if (j > 0)				Node[k][i][j]->LinkNodeXL(Node[k][i][j - 1]);
					//X-RIGHT
					if (j < nxNode - 1)		Node[k][i][j]->LinkNodeXR(Node[k][i][j + 1]);
					//Y-LEFT
					if (i > 0)				Node[k][i][j]->LinkNodeYL(Node[k][i - 1][j]);
					//Y-RIGHT
					if (i < nyNode - 1)		Node[k][i][j]->LinkNodeYR(Node[k][i + 1][j]);
					//Z-LEFT
					if (k > 0)				Node[k][i][j]->LinkNodeZL(Node[k - 1][i][j]);
					//Z-RIGHT
					if (k < nzNode - 1)		Node[k][i][j]->LinkNodeZR(Node[k + 1][i][j]);
				}
			}
		}
	}
}

void CStructure::LinkCXTable(CCXManager *pCXMan)
{
	int i, j, k;

	for (i = 0; i < nyNode; i++) {
		for (j = 0; j < nxNode; j++) {

			//Node[i][j]->CXs = new CCXTable*[nzNode];
			for (k = 0; k < nzNode; k++) {
				if (Node[k][i][j]) {

					Node[k][i][j]->SetCXTable(pCXMan->GetCXTable(Node[k][i][j]->GetCXTableID()));

					//Node[i][j]->CXs[k] = (pCXMan->GetCXTable(nodeTypeId[Node[i][j]->configId][k]));
				}
				//Node[i][j]->SetCXTable(pCXMan->GetCXTable(Node[i][j]->GetCXTableID()));
			//}
			}
		}
	}
}


int CStructure::CalculateNumberOfMeshes(real * MeshSize, int ** nMeshX, int * nMeshXY)
{
	int nx, ny, nz;
	int mnx, mnxy, mnxyz = 0;
	int i, j, k;
	int tot = 0;

	for (k = 0; k < nzNode; k++) {
		mnxy = 0;
		nz = (int)(zNodeSize[k] / MeshSize[ZDIR]);
		for (j = 0; j < nyNode; j++) {
			mnx = 0;
			ny = (int)(yNodeSize[j] / MeshSize[YDIR]);
			for (i = 0; i < nxNode; i++) {
				if (Node[k][j][i]) {
					nx = (int)(xNodeSize[i] / MeshSize[XDIR]);
					mnx += nx;
					Node[k][j][i]->SetMeshNumbers(nx, ny, nz);
				}
			}
			nMeshX[k][j] = mnx;
			mnxy += mnx*ny;
		}
		nMeshXY[k] = mnxy;
		mnxyz += mnxy*nz;
	}

	return mnxyz;
}
void CStructure::SetCoefficient(real * MeshSize,int nMeshAll)
{

	int i, j, k;
	int dir;
	for (k = 0; k < nzNode; k++) {
		for (j = 0; j < nyNode; j++) {
			for (i = 0; i < nxNode; i++) {
				if (Node[k][j][i])
				{
					for (dir = 0; dir < 3; dir++) {
						Node[k][j][i]->meshSize[dir] = MeshSize[dir];
					}
					Node[k][j][i]->SetCoefficient();
					
					
					Node[k][j][i]->Flux[0] = 1. / nMeshAll / 2;
					Node[k][j][i]->Flux[1] = 1. / nMeshAll / 2;

				}
			}
		}
	}
}

int CStructure::ExecutePowerIterationMethod(real *albedo, real keff,CVector &FluxOut, int ** nMeshX, int * nMeshXY)
{
	int i, j, k;
	//int ny = 0, nz = 0;
	//int loc[3] = { 0,0,0 };
	int numFlux = 0;
	int nIt;
	bool bConv = false;
	real maxErr;
	real prevKeff;

	real denom, norm;
	CVector preFlux;
	preFlux.SetDimension(FluxOut.GetDimension());

	CVector oneVector;
	oneVector.SetDimension(FluxOut.GetDimension());

	real sumSource;
	nIt = 0;

	real ktemp = 1.;

	do {
		preFlux.SetValues(FluxOut);

		prevKeff = keff;

		numFlux = 0;
		sumSource = 0;
		for (k = 0; k < nzNode; k++) {
			for (j = 0; j < nyNode; j++) {
				for (i = 0; i < nxNode; i++) {
					if (Node[k][j][i])
					{
 						/*Node[k][j][i]->Flux[0] = FluxOut.GetValue(numFlux*NUM_GRP);
						Node[k][j][i]->Flux[1] = FluxOut.GetValue(numFlux*NUM_GRP+1);*/
						Node[k][j][i]->Nem(albedo, keff);
						FluxOut.SetValue(numFlux*NUM_GRP, Node[k][j][i]->Flux[0]);
						FluxOut.SetValue(numFlux*NUM_GRP + 1, Node[k][j][i]->Flux[1]);

						numFlux++;

						sumSource += Node[k][j][i]->source[0] + Node[k][j][i]->source[1];
					}
				}
			}
		}

		denom = FluxOut.Multiply(preFlux);
		norm = FluxOut.Multiply(FluxOut);
	//	ktemp *= norm / denom;

	//	cout << ktemp << endl;
		//cout << sumSource << endl;
			
		keff = prevKeff*norm / denom;
		//keff = sumSource;
//		keff = FluxOut.Multiply(oneVector);

		maxErr = FluxOut.GetMaxRelDiff(preFlux);
		if (maxErr<1.e-6) bConv = true;

		for (k = 0; k < nzNode; k++) {
			for (j = 0; j < nyNode; j++) {
				for (i = 0; i < nxNode; i++) {
					if (Node[k][j][i])
					{
						Node[k][j][i]->Flux[0] /= keff;
						Node[k][j][i]->Flux[1] /= keff;
						numFlux++;
					}
				}
			}
		}
		nIt++;

		cout << " " << "Iteration : " << nIt << " " << " keff :" << "\t" << keff;
		cout << " Error = " << maxErr << endl;
	} while (!bConv);

	ofstream Flux11("Flux.txt");
	Flux11.precision(5);
	Flux11.setf(ios::scientific, ios::floatfield);
	for (k = 0; k < nzNode; k++) {
		Flux11 << "z: " << k + 1 << endl;
		for (j = 0; j < nyNode; j++) {
			for (i = 0; i < nxNode; i++) {
				if (Node[k][j][i])
				{
					Flux11 << setw(14) << Node[k][j][i]->Flux[0]<<'\t';
				}
			}
			Flux11 << endl;
		}
		Flux11 << endl<<endl;
	}

	ofstream Flux2("Flux2.txt");
	Flux2.precision(5);
	Flux2.setf(ios::scientific, ios::floatfield);
	for (k = 0; k < nzNode; k++) {
		Flux2 << "z: " << k + 1 << endl;
		for (j = 0; j < nyNode; j++) {
			for (i = 0; i < nxNode; i++) {
				if (Node[k][j][i])
				{
					Flux2 << setw(14) << Node[k][j][i]->Flux[1] << '\t';
				}
			}
			Flux2 << endl;
		}
		Flux2 << endl << endl;
	}



	Flux11.close();
	Flux2.close();


	return nIt;
}


//
//void CStructure::SetNetLossMatrix(CSparceMatrix * pMatL, real * MeshSize, int ** nMeshX, int * nMeshXY, int nMeshAll, real * albedo)
//{
//	int i, j,k ;
//	int ny = 0, nz = 0;
//	int loc[3] = { 0,0,0 };
//
//	for (k = 0; k < nzNode; k++) {
//		for (j = 0; j < nyNode; j++) {
//			for (i = 0; i < nxNode; i++) {
//				if (Node[k][j][i]) {
//					Node[k][j][i]->SetNetLossMatrix(pMatL, MeshSize, loc, k, j, nMeshX, nMeshXY, nMeshAll, albedo);
//					loc[0] += Node[k][j][i]->GetNumberOfMesh(XDIR);
//				}
//				if (i == 0) {
//					ny = Node[k][j][i]->GetNumberOfMesh(YDIR);
//					nz = Node[k][j][i]->GetNumberOfMesh(ZDIR);
//				}
//			}
//			loc[0] = 0;
//			loc[1] += nMeshX[k][j] * ny;
//		}
//		loc[1] = 0;
//		loc[2] += nMeshXY[k] * nz;
//	}
//}
//
//void CStructure::SetFissProdMatrix(CSparceMatrix * pMatF, real * MeshSize, int ** nMeshX, int * nMeshXY, int nMeshAll)
//{
//	int i, j,k ;
//	int ny = 0, nz = 0;
//	int loc[3] = { 0,0,0 };
//
//	for (k = 0; k < nzNode; k++) {
//		for (j = 0; j < nyNode; j++) {
//			for (i = 0; i < nxNode; i++) {
//				if (Node[k][j][i]) {
//					Node[k][j][i]->SetFissProdMatrix(pMatF, MeshSize, loc, k, j, nMeshX, nMeshXY, nMeshAll);
//					loc[0] += Node[k][j][i]->GetNumberOfMesh(XDIR);
//				}
//				if (i == 0) {
//					ny = Node[k][j][i]->GetNumberOfMesh(YDIR);
//					nz = Node[k][j][i]->GetNumberOfMesh(ZDIR);
//				}
//			}
//			loc[0] = 0;
//			loc[1] += nMeshX[k][j] * ny;
//		}
//		loc[1] = 0;
//		loc[2] += nMeshXY[k] * nz;
//	}
//}
//
//void CStructure::SetXiMatrix(CSparceMatrix & Xi, real * MeshSize, int ** nMeshX, int * nMeshXY, int nMeshTot)
//{
//	int i, j, k;
//	int ny = 0, nz = 0;
//	int loc[3] = { 0,0,0 };
//
//	for (k = 0; k < nzNode; k++) {
//		for (j = 0; j < nyNode; j++) {
//			for (i = 0; i < nxNode; i++) {
//				if (Node[k][j][i]) {
//					Node[k][j][i]->SetXiMatrix(Xi, MeshSize, loc, nMeshX[k][j], nMeshXY[k], nMeshTot);
//					loc[0] += Node[k][j][i]->GetNumberOfMesh(XDIR);
//				}
//				if (i == 0) {
//					ny = Node[k][j][i]->GetNumberOfMesh(YDIR);
//					nz = Node[k][j][i]->GetNumberOfMesh(ZDIR);
//				}
//			}
//			loc[0] = 0;
//			loc[1] += nMeshX[k][j] * ny;
//		}
//		loc[1] = 0;
//		loc[2] += nMeshXY[k] * nz;
//	}
//}
