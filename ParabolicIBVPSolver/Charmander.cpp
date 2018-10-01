#include "tridiagonal_matrix.h"
#include "twoPointBVP.h"
#include "TwoPointBVPAppr.h"
#include "ParabolicIBVPSolver.h"
#include <iostream>
#include <fstream>

//-----------------------------------------------
// declare any needed constants
//-----------------------------------------------

double pi = acos(-1.0);

int numspace = 50;

int numtime =50;

double endtime = 10.0;

//--------------------------------------------------
//set the diffusion coeffcient and if present 
// the reaction, forcing function and true solution
//--------------------------------------------------

double diffusioncoeff(vector<double> &x)
{
	return 1.0;
}

double forcecoeff(vector <double> &par)
{
	double x = par[0];
	double t = par[1];
	return exp(-t)*(0.5 *(1.0 - t)*(1.0 - x * x) + t);
}

double reactioncoeff(vector<double> &par)
{
	return 0.0;
}

double dudr(vector<double> &par)
{
	return -1.0;
}

double leftbdryvalue(double &t)
{
	return exp(-t)*cos(-1.0) - 1.0;
}

double rightbdryvalue(double &t)
{
	return exp(-t)*cos(1.0) + 1.0;
}

double initialcondition(double x)
{
	return x + cos(x);
}


double trusol(vector<double> &vars)
//vars[0] is space, vars[1] is time
{
    double x = vars[0];
    double t = vars[1];
    return exp(-t)*( cos(x) +0.5*t*(1 - x*x) ) + x;
}


int main()
{
	//---------------------------------------------------------------------
	//set up the two point bvp
	//---------------------------------------------------------------------
	double * dom = new double[2];
	dom[0] = -1.0;
	dom[1] = 1.0;

	TwoPointBVP *prob = new TwoPointBVP(dom, diffusioncoeff);

	double gammaa = 0.0;
	prob->set_left_bdry(true, gammaa, leftbdryvalue);

	double gammab = 0.0;

	prob->set_right_bdry(true, gammab, rightbdryvalue);

	prob->set_forcing_function(forcecoeff);
    
    prob->set_true_solution(trusol);

	//---------------------------------------------------------------------
	//Create the spatial the discritaization
	//---------------------------------------------------------------------

	//Declare the number of spatial subintervals and set up each subinterval
	int numsubintervals = numspace;

	double * subintervals = new double[numsubintervals];

	for (int i = 0; i < numsubintervals; i++)
	{
		subintervals[i] = (dom[1] - dom[0]) / numsubintervals;
	}

	//Create BVP Approximation and begin setting up spatial semi-discretization
	TwoPointBVPAppr * apprbvp = new TwoPointBVPAppr(numsubintervals,
		subintervals, prob);

	vector<double> xcoord = apprbvp->get_xcoord();

	//begin setting up the time full discretization by partioning time 
	int numtimesteps = numtime;
	double finaltime = endtime;
	double timestep = finaltime / numtimesteps;
	vector<double> timelevel(numtimesteps + 1);
	timelevel[0] = 0.0;

	// for loop to fill time-steps vector
    for (int i = 1; i <= numtimesteps; i++)
	{
		timelevel[i] = timelevel[i - 1] + timestep;
	}
    
    prob->AssembleBdrys(timelevel[0]);
    
    
    //display some info about the two point bvp
    //---------------------------------------------------------------------
    prob->display_info_TwoPointBVP();

    
	//Crete the Finite Diff method object
	ParabolicIBVPSolver *fdmethod;

	//Specify which method is used
	fdmethod = new BackwardEuler(apprbvp);

	//create the time interval and intial guess
	vector<double> timeinterval(2);
	vector<double> Ul(numsubintervals + 1);
	for (int i = 0; i <= numsubintervals; i++)
	{
		Ul[i] = initialcondition(xcoord[i]);
	}
    
    //Create a vector of vectors to store the solution.
    vector <vector<double> > Solns;
    
	//Create Ur vector(solution at  right time boundary) and iteratively
	//solve for it and store it in solution at each time interval
    //from t_0 to t_F
	
    vector<double> Ur;
    
    //Store the intial condition in the first vector
    Solns.push_back(Ul);
    
	for (int i = 1; i <= numtimesteps; i++)
	{
		timeinterval[0] = timelevel[i - 1];
		timeinterval[1] = timelevel[i];
		Ur = fdmethod->one_step_march(timeinterval, Ul);
        Solns.push_back(Ur);
        Ul = Ur;
	}

    //Create a true solution matrix
    vector< vector<double> > trueSolution(numtimesteps+1,
		vector<double>(numsubintervals+1) );
   
    vector<double> var(2);
    
    //fill the true solution matrix.
    for (int i = 0; i<numtimesteps+1; i++) {
        for (int j = 0; j<numsubintervals+1; j++){
            var[1] = timelevel[i];
            var[0] = xcoord[j];
          trueSolution[i][j] = prob->eval_true_solution(var);
        }
    }
    
    
    
    //output the approximate solution to a file
    ofstream outputFile;
    outputFile.open("approxSoln.txt");
    
    for (int i = 0; i<numtimesteps+1; i++){
        for (int j = 0; j<numsubintervals+1; j++){
            outputFile << Solns[i][j] << " ";
        }
        outputFile << endl;
    }
    outputFile.close();
    
    
    //output the true solution to a file

    outputFile.open("trueSoln.txt");
    
    for (int i = 0; i<numtimesteps+1; i++){
        for (int j = 0; j<numsubintervals+1; j++){
            outputFile << trueSolution[i][j] << " ";
        }
        outputFile << endl;
    }
    outputFile.close();
    
    //find the error at each point and ouput it to a file AprroxmationError
    vector< vector<double> >approximationError(numtimesteps+1,
		vector<double>(numsubintervals +1) );
   
    outputFile.open("approximationError.txt");
    
    for (int i = 0; i<numtimesteps+1; i++){
        for (int j = 0; j<numsubintervals+1; j++){
            outputFile << trueSolution[i][j] - Solns[i][j] << " ";
            approximationError[i][j] = trueSolution[i][j] - Solns[i][j];
        }
        outputFile << endl;
    }
    outputFile.close();
    
    // create a second derivative matrix and multiply it with
    // the solution vector.
    tridiagonal_matrix *scndDerv;
    scndDerv = new tridiagonal_matrix(numsubintervals +1);
    
    // fill first row
    scndDerv->add_to_diagonal_entry(0, 0);
    scndDerv->add_to_upper_diagonal_entry(0, 0);
    
    //fill interior
    for (int i = 1; i<numsubintervals; i++) {
        scndDerv->add_to_upper_diagonal_entry(i, 1.0);
        scndDerv->add_to_diagonal_entry(i, -2.0);
        scndDerv->add_to_lower_diagonal_entry(i-1, 1.0);
    }
    
    //fill last row
    //scndDerv->add_to_lower_diagonal_entry(numsubintervals-1, 0.0);
    //scndDerv->add_to_diagonal_entry(numsubintervals, 0.0);
    
    //Create matrix to store second derv of solution U.
    vector< vector<double> >Derv2;
    
    outputFile.open("2DervSurf.txt");
    
    //for loop to find d_xx of U for t=0 to t=final time
    // and then ouptut it to a file
    for (int i=0; i<numtimesteps+1; i++) {
        double h = (dom[1] - dom[0])/numsubintervals;
        vector<double> derv(numsubintervals+1);
        derv=scndDerv->Mult(Solns[i]);
        derv[0] = (Solns[i][0] - 2*(Solns[i][1]) + Solns[i][2])/h;
        derv[numsubintervals] =(Solns[i][numsubintervals-2]
                - 2*(Solns[i][numsubintervals-1])
                + Solns[i][numsubintervals])/h;
        Derv2.push_back(derv);
        for (int j=0; j < numsubintervals+1; j++) {
            outputFile << derv[j] <<" ";
        }
        outputFile << endl;
    }
    outputFile.close();
    
    
	//Clean up memory
	delete apprbvp;
	delete fdmethod;
	delete[]subintervals;
	delete prob;
	delete[]dom;
    delete scndDerv;

	return 0;
}
