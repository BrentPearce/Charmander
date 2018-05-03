#include "tridiagonal_matrix.h"
#include "twopointbvp.h"
#include "twopointbvpappr.h"
#include "ParabolicIBVPSolver.h"


ParabolicIBVPSolver::ParabolicIBVPSolver(TwoPointBVPAppr * prob)
{
	apprbvp = prob;
	diffusionmat = apprbvp->calcDiffusion(); // memory is created.
}

ParabolicIBVPSolver::~ParabolicIBVPSolver()
{
	delete diffusionmat;
}

//----------------------------------------------------------------------------

ForwardEuler::ForwardEuler(TwoPointBVPAppr * prob) : ParabolicIBVPSolver(prob)
{
	cout << " A Forward Euler Scheme is utilized\n";
}

vector<double> ForwardEuler::one_step_march(vector<double> &timeinterval,
	vector<double> &Ul)
{
	vector<double> Ur;
	// implement here according the derivation we did in class (see class notes).

	return Ur;
}

//----------------------------------------------------------------------------

BackwardEuler::BackwardEuler(TwoPointBVPAppr * prob) : ParabolicIBVPSolver(prob)
{
	cout << " A Backward Euler Scheme is ultilized\n";

}

double find_l2_norm(vector<double> const x)
{
	int num_entries = x.size();
	double sum_o_squares = 0;
	for (int i = 0; i < num_entries; i++)
	{
		sum_o_squares = sum_o_squares + x[i] * x[i];
	}
	return sqrt(sum_o_squares);
}

vector<double> BackwardEuler::one_step_march(vector<double> &timeinterval, 
											 vector<double> &Ul)
{
	//create some needed vars
	int numsubintervals = apprbvp->get_numsubintervals();
	int iteration_counter = 0;
	int max_iterations = 100;
	double Toler = 1e-20;
	double norm =1e10;
	tridiagonal_matrix *Gp, *A;
	vector<double> R(numsubintervals + 1, 0.0);
	vector<double>Rp(numsubintervals + 1, 0.0);
	vector<double>F;
	double deltaTime = timeinterval[1] - timeinterval[0];
	A = new tridiagonal_matrix(numsubintervals + 1);
	// Calculate the tridiagonal matrix coming from diffusion component.
	A = apprbvp->calcDiffusion();
	vector<double> U_l = Ul;
	vector<double> Ur = Ul;

	vector<double> lumpMassMatrix = apprbvp->calcLumpedMass();

	vector<double> lumpMatrixInvers(numsubintervals + 1);
	


	while (iteration_counter <= max_iterations && norm>Toler)
	{
		Gp = new tridiagonal_matrix(A);
		


	}


	return Ur;
}

//----------------------------------------------------------------------------

ImplicitMidPoint::ImplicitMidPoint(TwoPointBVPAppr * prob) : ParabolicIBVPSolver(prob)
{
	cout << " An Implicit MidPoint Scheme is utilized\n";
}

vector<double> ImplicitMidPoint::one_step_march(vector<double> &timeinterval,
	vector<double> &Ul)
{
	vector<double> Ur;

	for (int i = 0; i <= numtimeintervals; i++)
	{
		Ur[i] = Ul[i] - ((timeinterval[1] - timeinterval[0]) * MInv[i] * G((1 / 2)*(timeinterval[1] + timeinterval[0]),
			(1 / 2)*(Ul[i] + Ur[i])));
	}

	return Ur;
}

//----------------------------------------------------------------------------

ExplicitMidPoint::ExplicitMidPoint(TwoPointBVPAppr * prob) : ParabolicIBVPSolver(prob)
{
	cout << " An Explicit MidPoint Scheme is utilized\n";
}

vector<double> ExplicitMidPoint::one_step_march(vector<double> &timeinterval,
	vector<double> &Ul)
{
	vector<double> Ur;
	vector<double> Um;
	double tm = (1 / 2) * (timeinterval[1] + timeinterval[0]);

	for (int i = 0; i <= numtimesteps; i++)
	{
		Um[i] = Ul[i] - ((tm - timeinterval[0]) * MInv[i] * G(timeinterval[0], Ul[i]));
		Ur[i] = Ul[i] - ((timeinterval[1] - timeinterval[0]) * G((1 / 2)*(timeinterval[1] - timeinterval[0]), Um[i]));
	}

	return Ur;
}