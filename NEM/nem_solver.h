#ifndef _nem_SOLVER_H_
#define _nem_SOLVER_H_

#include "define.h"
#include "cxman.h"
#include "node.h"

#include "vector.h"

/**************************************************************************************************
 * nem Solver
**************************************************************************************************/
class CnemSolver {
//-------------------------------------------------------------------------------------------------
// Objects for the core calculations
//-------------------------------------------------------------------------------------------------
public:
    CCXManager CXMan; //CX Manager object
    CStructure Str; //Structure object

	CSparceMatrix MatL; //Net-Loss Matrix
	CSparceMatrix MatF; //Fission Matrix
	CSparceMatrix Xi;
//-------------------------------------------------------------------------------------------------
// Calculation Options
//-------------------------------------------------------------------------------------------------
private:
    int nDim;
    real meshSize[3]; //in the order of XDIR, YDIR, ZDIR
	real albedo[6];
    real ke; //expected keff for the Wieland method
    real omega; //extrapolation parameter for the SOR method

    int nHigherOrder; //the maximum order of eigenvectors to calculate
	int method;

//-------------------------------------------------------------------------------------------------
// Output
//-------------------------------------------------------------------------------------------------
private:
    real keff;
    CVector Flux;

public:
    CnemSolver(void); //Constructor
    ~CnemSolver(void); //Destructor

//-------------------------------------------------------------------------------------------------
// Methods
//-------------------------------------------------------------------------------------------------
public:
    void ReadInput(char *filename); //Read Input file
        void ReadOptions(istream &ins); //Read "Options" section in a input file

	void AnalyzeCore(void); //Analyze the core problem by the nem method
    int ExecutePowerIterationMethod(CSparceMatrix &MatL, CSparceMatrix &MatF, CVector &Flux, real &keff);

};

#endif