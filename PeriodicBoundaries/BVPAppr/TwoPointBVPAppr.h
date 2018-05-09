#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>


using namespace std;

class TwoPointBVPAppr
{
protected:

	int numsubintervals;
	//the number of intervals in the domain

	const double *steplenghts;
	const TwoPointBVP *theproblem;
	vector <double> xcoord;
	vector <double> midcoord;
	vector <double> Deltax;
	bool guess_seed_is_present;

	double AssembleDiffusion(tridiagonal_matrix * tmat);

	void AssembleReaction(vector<double> &W, vector<double> &RW, 
							vector<double> &RPW);

	vector<double> AssembleForce();

	double(*intialGuessSeed)(vector<double> &);

public:

	TwoPointBVPAppr(int N, const double *subintervallengths, const TwoPointBVP *prob);

	vector<double> get_xcoord();

	int get_numsubintervals();

	void set_intial_guess_seed(double(*guessSeed)(vector<double> &));

	double eval_intial_guess_seed(vector<double> &par) const;

	vector<double> Solve(int max_num_iter, double TOL);

	double find_max_error(int max_iters, double TOL);

	tridiagonal_matrix * calcDiffusion();



	vector<double> calcLumpedMass();


	~TwoPointBVPAppr();
};
