#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <cmath>

using namespace std;

// Base/parent class for solving Parabolic IBVP
class ParabolicIBVPSolver {

protected:
	TwoPointBVPAppr * apprbvp;
	tridiagonal_matrix * diffusionmat;

public:

	// class constructor
	ParabolicIBVPSolver(TwoPointBVPAppr * prob);

	/* This is a pure virtual function, i.e., the parent does not have
	an implementation. Its actual implementation is done in the child
	class and is specific to a child.
	*/

	virtual vector<double> one_step_march(vector<double> &timeinterval,
		vector<double> &Ul) = 0;

	// class destructor
	virtual ~ParabolicIBVPSolver();
};

//---------------------------------------------------------------------------

class ForwardEuler : public ParabolicIBVPSolver {

private:

public:

	ForwardEuler(TwoPointBVPAppr * prob);

	virtual vector<double> one_step_march(vector<double> &timeinterval,
		vector<double> &Ul);
};

//---------------------------------------------------------------------------

class BackwardEuler : public ParabolicIBVPSolver {

private:

public:

	BackwardEuler(TwoPointBVPAppr * prob);

	virtual vector<double> one_step_march(vector<double> &timeinterval,
		vector<double> &Ul);
};

//---------------------------------------------------------------------------

class ImplicitMidPoint : public ParabolicIBVPSolver {

private:

public:

	ImplicitMidPoint(TwoPointBVPAppr * prob);

	virtual vector<double> one_step_march(vector<double> &timeinterval,
		vector<double> &Ul);
};

//---------------------------------------------------------------------------

class ExplicitMidPoint : public ParabolicIBVPSolver {

private:

public:

	ExplicitMidPoint(TwoPointBVPAppr * prob);

	virtual vector<double> one_step_march(vector<double> &timeinterval,
		vector<double> &Ul);
};
