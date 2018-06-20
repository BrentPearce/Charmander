#include "tridiagonal_matrix.h"
#include "twopointbvp.h"
#include "twopointbvpappr.h"
#include "ParabolicIBVPSolver.h"
#include <iostream>
#include <fstream>

//-----------------------------------------------
// declare any needed constants
//-----------------------------------------------

double pi = acos(-1.0);

double lambda = 2.0;

//double theta = 2.3575510539e+00;

double theta = 8.5071995707e+00;

int numspace = 5000;

int numtime = 100;

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

//double reactioncoeff(vector<double> &par)
//{
//	return 0.0;
//}

//double dudr(vector<double> &par)
//{
//	return -1.0;
//}

double leftbdryvalue(double &t)
{
	
	return exp(-t)*cos(1.0) - 1.0;
}

double rightbdryvalue(double &t)
{
	return exp(-t)*cos(1.0) + 1.0;
}

double initialcondition(double x)
{
	return x + cos(x);
}

bool true_sol_is_present = true;

double trusol(double x, double t)
{
	return exp(-t)*(cos(x) + 0.5*t*(1.0 - x * x)) + x;
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

	
	//---------------------------------------------------------------------
	//display some info about the two point bvp
	//---------------------------------------------------------------------
	prob->display_info_TwoPointBVP();

	//---------------------------------------------------------------------
	//Run the aproximation
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
	
	for (int i = 1; i <= numtimesteps; i++)
	{
		timelevel[i] = timelevel[i - 1] + timestep;
	}
	//Crete the FinDiff method object
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
	//Create Ur vector(solution at  right time boundary) and iteratively
	//solve for it at each time interval from t_0 to t_F
	vector<double> Ur;
	for (int i = 1; i <= numtimesteps; i++)
	{
		timeinterval[0] = timelevel[i - 1];
		timeinterval[1] = timelevel[i];
		Ur = fdmethod->one_step_march(timeinterval, Ul);
		Ul = Ur;
	}

	//Clean up memory
	delete apprbvp;
	delete fdmethod;
	delete[]subintervals;
	delete prob;
	delete[]dom;


	return 0;
}
