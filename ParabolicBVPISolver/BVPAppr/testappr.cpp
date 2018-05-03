#include "tridiagonal_matrix.h"
#include "twopointbvp.h"
#include "twopointbvpappr.h"
#include <iostream>
#include <fstream>

//-----------------------------------------------
// declare any needed constants
//-----------------------------------------------

double pi = acos(-1.0);

double lambda = 2.0;

//double theta = 2.3575510539e+00;

double theta = 8.5071995707e+00;

//double phi = 0.0;

int numsub = 100;

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

double leftbdryvalue(vector<double> &par)
{
	double t = par[0];
	return exp(-t)*cos(1.0) - 1.0;
}

double rightbdryvalue(vector<double> &par)
{
	double t = par[0];
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
	//------------------------------------------------
	//set up the two point bvp
	//------------------------------------------------
	double * dom = new double[2];
	dom[0] = -1.0;
	dom[1] = 1.0;

	TwoPointBVP *prob = new TwoPointBVP(dom, diffusioncoeff);

	double gammaa = 0.0;
	prob->set_left_bdry(true, gammaa, leftbdryvalue);

	double *rbval = new double[2];
	rbval[0] = 0.0; //this is gamma_l
	rbval[1] = 0.0;//this is g_l
	prob->set_right_bdry(true, rbval);

	
	prob->set_forcing_function(forcecoeff);

	prob->set_true_solution(truesol);


	//-----------------------------------------------------
	//display some info about the two point bvp
	//-----------------------------------------------------
	prob->display_info_TwoPointBVP();



	return 0;
}
