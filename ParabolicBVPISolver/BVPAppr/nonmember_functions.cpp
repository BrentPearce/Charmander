#include"nonmember_functions.h"

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