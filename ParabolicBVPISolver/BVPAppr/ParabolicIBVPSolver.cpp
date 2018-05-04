#include "tridiagonal_matrix.h"
#include "twopointbvp.h"
#include "twopointbvpappr.h"
#include "ParabolicIBVPSolver.h"

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
	std::ofstream ofs;
	ofs.open("problem_info.txt", std::ofstream::out | std::ofstream::app);
	ofs << " A Backward Euler Scheme is ultilized\n";
	ofs.close();
}

vector<double> BackwardEuler::one_step_march(vector<double> &timeinterval, 
											 vector<double> &Ul)
{
	//create some needed vars
	int numsubintervals = apprbvp->get_numsubintervals();
	int iteration_counter = 0;
	int max_iterations = 100;
	double Toler = 1e-20;
	double norm = 1e10;
	tridiagonal_matrix *Gp, *A;
	vector<double> R(numsubintervals + 1, 0.0);
	vector<double>Rp(numsubintervals + 1, 0.0);
	vector<double>F;
	double deltaTime = timeinterval[1] - timeinterval[0];
	A = new tridiagonal_matrix(numsubintervals + 1);
	// Calculate the tridiagonal matrix coming from diffusion component.
	A = apprbvp->calcDiffusion();
	
	//create the left time solution and intial right time solution
	vector<double> U_l = Ul;
	vector<double> Ur = Ul;

	

	// fill the force vector
	F = apprbvp->calcForce(timeinterval[1]);

	vector<double> lumpMassMatrix = apprbvp->calcLumpedMass();

	
	vector<double> G(numsubintervals +1);

	//Newton Iteration solving Gp*h = G 
	while (iteration_counter <= max_iterations && norm>Toler)
	{
		Gp = new tridiagonal_matrix(A);
		
		vector<double> AUr = A->Mult(Ur);

		// if there is a reaction calculate r(x,u) and pd(r(x,u),u)
		// (if there isn't then fill R and RP all zero)
		apprbvp->calcReaction(Ur, R, Rp);

		for (int i = 0; i < numsubintervals + 1; i++)
		{
			Gp->add_to_diagonal_entry(i, Rp[i] +
				(1 / deltaTime)*lumpMassMatrix[i]);
		}
		for (int i = 0; i < numsubintervals + 1; i++)
		{
			G[i] = -1*((1 / deltaTime)*lumpMassMatrix[i]*(Ur[i] - Ul[i])) -
					(AUr[i] + R[i] -F[i]);
		}
		
		for (int i = 0; i < numsubintervals + 1; i++)
		{
			Gp->add_to_diagonal_entry(i, (Rp[i] +
				(1 / deltaTime)*lumpMassMatrix[i]));
		}
		
		//transform Gp so we can call solve_linear system
		Gp->transform();

		//solve for delta needed to update 
		vector<double> delta(numsubintervals + 1);
		delta = Gp->solve_linear_system(G);

		//Update Ur
		for (int i = 0; i <= numsubintervals; i++)
		{
			Ur[i] = Ur[i] + delta[i];
		}

		//delete the Tridiagonal Matrix Gp associated with the iteration
		delete Gp;

		//find the norm of delta to see if iterations continue
		norm = find_l2_norm(delta);
		
		//update the iteration counter.
		iteration_counter++;

	}

	if (iteration_counter == max_iterations)
	{
		std::ofstream ofs;
		ofs.open("problem_info.txt", std::ofstream::out | std::ofstream::app);

		ofs << " Convergence not reached within max number of iterations:  "
			<< max_iterations << endl;

		ofs.close();
	}
	

	return Ur;
}

//-------------------------------------------------------------------------

ImplicitMidPoint::ImplicitMidPoint(TwoPointBVPAppr*prob):ParabolicIBVPSolver(prob)
{
	cout << " An Implicit MidPoint Scheme is utilized\n";
}

vector<double> ImplicitMidPoint::one_step_march(vector<double> &timeinterval,
	vector<double> &Ul)
{
	vector<double> Ur;



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
	
	return Ur;
}