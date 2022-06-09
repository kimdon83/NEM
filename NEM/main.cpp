#include "nem_solver.h"

#include "define.h"
int main(void)
{
    CnemSolver nemSol;

	//time start
	clock_t start_point, end_point;
	start_point = clock();
	double caltime;

    //Read an input
//	nemSol.ReadInput("input_1D_2G_hom_test.txt");      //	0.884768

//	nemSol.ReadInput("input_IAEA1.txt"); 			//  1.02911
//	nemSol.ReadInput("input_IAEA1-1.txt"); 			//  1.02912
	nemSol.ReadInput("input_IAEA2.txt"); 			//  1.02908

    //Solve the core analysis problem
	nemSol.AnalyzeCore();

	end_point = clock();
	caltime = (end_point - start_point)/(double)( CLOCKS_PER_SEC);

	cout << caltime << "(sec)" << endl;

	system("pause");
    return 0;
	string abc;
}
