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

vector<double> BackwardEuler::one_step_march(vector<double> &timeinterval, vector<double> &Ul)
{
	tridiagonal_matrix *A;
	A = diffusionmat;
	vector<double> F;

	for (int i = 0; i <= numtimeintervals; i++)
	{
		Ur[i] = Ul[i] - (timeinterval[1] - timeinterval[0]) * MInv[i] * G[i];
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