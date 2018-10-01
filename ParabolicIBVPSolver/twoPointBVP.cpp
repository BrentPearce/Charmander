#include "twoPointBVP.h"


// A C++ class implementation for two point BVPs and some other stuff

TwoPointBVP::TwoPointBVP(double *dom, double(*dFunc) (vector<double> &))
{
	domain = dom;
	diffusion = dFunc;
	reactionIsPresent = false;
	forcingFunctIsPresent = false;
}

void TwoPointBVP::set_left_bdry(bool _leftIsDirichlet, double *val)
{
	leftBdryIsDirichlet = _leftIsDirichlet;
	leftBdryValues = val;
}

void TwoPointBVP::set_right_bdry(bool _rightIsDirichlet, double *val)
{
	rightBdryIsDirichlet = _rightIsDirichlet;
	rightBdryValues = val;
}

void TwoPointBVP::set_reaction(double(*rFunctOne) (vector<double> &),
	double(*rFunctTwo) (vector<double> &))
{
	reaction = rFunctOne;
	partialreactionpartialu = rFunctTwo;
	reactionIsPresent = true;
}

void TwoPointBVP::set_forcing_function(double(*fFunct) (vector<double> &))
{
	forcingFunct = fFunct;
	forcingFunctIsPresent = true;
}

void TwoPointBVP::set_true_solution(double(*TrueSol)(vector<double>&))
{
	trueSolu = TrueSol;
	trueSolIsPresent = true;
}

double * TwoPointBVP::get_domain() const
{
	return domain;
}

bool TwoPointBVP::left_bdry_is_Dirichlet() const
{
	return leftBdryIsDirichlet;
}

bool TwoPointBVP::right_bdry_is_Dirichlet() const
{
	return rightBdryIsDirichlet;
}

double * TwoPointBVP::get_left_bdry_values() const
{
	return leftBdryValues;
}

double * TwoPointBVP::get_right_bdry_values() const
{
	return rightBdryValues;
}

bool TwoPointBVP::reaction_is_present() const
{
	return reactionIsPresent;
}

bool TwoPointBVP::forcing_fucntion_is_present() const
{
	return forcingFunctIsPresent;
}

bool TwoPointBVP::true_solution_is_present() const
{
	return trueSolIsPresent;
}

double TwoPointBVP::eval_diffusion(vector<double> &x) const
{
	return diffusion(x);
}

vector<double> TwoPointBVP::eval_reaction(vector<double> &par) const
{
	vector<double> val(2);
	val[0] = reaction(par);
	val[1] = partialreactionpartialu(par);
	return val;
}

double TwoPointBVP::eval_forcing_function(vector<double>& x) const
{
	return forcingFunct(x);
}

double TwoPointBVP::eval_true_solution(vector<double> &x) const
{
	return trueSolu(x);
}

void TwoPointBVP::display_info_TwoPointBVP() const
{
	ofstream fileout;
	fileout.open("problem_info.txt");
	fileout << "************************************************************************ \n";
	fileout << " Some info regarding the two" <<
		" point BVP problem and approximation: \n";


	fileout << " Domain is:(" << domain[0] << "," << domain[1] << ")"
		<< endl;

	if (leftBdryIsDirichlet)
	{
		fileout << " Initial left boundary is Dirichlet with value: "
			<< leftBdryValues[1] << endl;
	}
	else if (leftBdryValues[0] == 0)
	{
		fileout << " Initial left Boundary is Neumann with g_0: "
			<< leftBdryValues[1] << endl;
	}
	else
	{
		fileout << " Initial left Boundary is Robin with gamma_0: "
			<< leftBdryValues[0]
			<< " and g_0: " << leftBdryValues[1] << endl;
	}

	if (rightBdryIsDirichlet)
	{
		fileout << " Initial right boundary is Dirichlet with value: "
			<< rightBdryValues[1] << endl;
	}
	else if (rightBdryValues[0] == 0)
	{
		fileout << " Initial right Boundary is Neumann with g_L: "
			<< rightBdryValues[1] << endl;
	}
	else
	{
		fileout << " Initial right Boundary is Robin with gamma_L: "
			<< rightBdryValues[0]
			<< " and g_L: " << rightBdryValues[1] << endl;
	}

    
	if (reactionIsPresent)
	{
		fileout << " Reaction is Present. \n";
	}
	else
	{
		fileout << " No reaction is present. \n";
	}

	if (forcingFunctIsPresent)
	{
		fileout << " Forcing Function is Present. \n";
	}
	else
	{
		fileout << " Forcing function is not present. \n";
	}

	fileout << "******************************************************* \n";
	fileout.close();
	return;
}



//---------------For PDEs--------------------------------------------------

//-------------------------------------------------------------------------

void TwoPointBVP::set_left_bdry(bool _leftIsDirichlet, double gamma_l,
	double(*g_a)(double &))
{
	leftBdryIsDirichlet = _leftIsDirichlet;
	ga = g_a;
	gamma_a = gamma_l;
}

void TwoPointBVP::set_right_bdry(bool _rightIsDirichlet, double gamma_r,
	double(*g_b)(double &))
{
	rightBdryIsDirichlet = _rightIsDirichlet;
	gb = g_b;
	gamma_b = gamma_r;
}

double TwoPointBVP::calcLeftBdry(double t) const
{
	return ga(t);
}

double TwoPointBVP::calcRightBdry(double t) const
{
	return gb(t);
}

void TwoPointBVP::AssembleBdrys(double &t)
{
	double *val1 = new double[2];
	val1[0] = gamma_a;
	val1[1] = calcLeftBdry(t);
	leftBdryValues = val1;
    
    double *val2 = new double[2];
    val2[0] = gamma_b;
	val2[1] = calcRightBdry(t);
	rightBdryValues = val2;
}



TwoPointBVP::~TwoPointBVP()
{
	;
}
