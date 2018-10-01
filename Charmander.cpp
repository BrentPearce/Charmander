#include "tridiagonal_matrix.h"
#include "twoPointBVP.h"
#include "TwoPointBVPAppr.h"
#include <iostream>
#include <fstream>

//-----------------------------------------------
// declare any needed constants
//-----------------------------------------------

double pi = acos(-1.0);

int numsub = 100;

int const maxIters = 100;
double const Toler = 1.0e-9;

//--------------------------------------------------
//set the diffusion coeffcient and if present 
// the reaction, forcing function and true solution
//--------------------------------------------------

double diffusioncoeff(vector<double> &x)
{
   return 1.0;
}

double forcecoeff(vector <double> &x)
{
    return sin(4*pi*x[0]);
}

double reactioncoeff(vector<double> &par)
{
    return 6.0*( cbrt(par[1]) );
}

double dudr(vector<double> &par)
{
    return 2.0 / (pow(par[1], 2.0/3.0) );
}

double truesol(vector<double> & x)
{
    return (1.0/16*pi*pi)*(sin(4*pi*x[0]));
}

double seed(vector<double> & par)
{
    return 0.0;
}

int main()
{
    //------------------------------------------------
    //set up the two point bvp
    //------------------------------------------------
    double * dom = new double[2];
    dom[0] = -2.0;
    dom[1] = 2.0;
    
    TwoPointBVP *prob = new TwoPointBVP(dom, diffusioncoeff);
    
    //prob->set_left_periodic_bdry();
    
    //prob->set_right_periodic_bdry();
    
    
    double *lbval = new double[2];
    lbval[0] = 0.0; // this is gamma_l
    lbval[1] = 0.0; // this is g_l
    prob->set_left_bdry(true, lbval);
    
    
    double *rbval = new double[2];
    rbval[0] = 0.0; //this is gamma_r
    rbval[1] = 0.0;//this is g_r
    prob->set_right_bdry(true, rbval);
    
    
    //prob->set_reaction(reactioncoeff,dudr);
    
    prob->set_forcing_function(forcecoeff);
    
    prob->set_true_solution(truesol);
    
	//-----------------------------------------------------
	//display some info about the two point bvp
	//-----------------------------------------------------
	prob->display_info_TwoPointBVP();


	//-----------------------------------------------------
	//solve for the approximate solution
	//(comment out this section if only finding error
	//higlight then ctrl + k, ctrl + c to comment out
	// highlight then ctrl+k, ctrl+u to uncomment)
	//-----------------------------------------------------
	
	{
		// create the 2 pt bvp approximation
		int numsubintervals = numsub;
		double * subintervals = new double[numsubintervals];
		for (int i = 0; i < numsubintervals; i++)
		{
			subintervals[i] = (dom[1] - dom[0]) / numsubintervals;
		}

		TwoPointBVPAppr *method = new TwoPointBVPAppr(numsubintervals,
			subintervals, prob);

		method->set_intial_guess_seed(seed);

		vector<double> sol = method->Solve(maxIters,Toler);
		vector<double> xcoord = method->get_xcoord();

        //Output the xcoords
        ofstream fileout;
        fileout.open("/Users/bpearce/octav/4/xcoords.txt");
        for (int i = 0; i < numsubintervals + 1; i++)
            fileout << xcoord[i] << endl;
        fileout.close();

        //Output the approximate Solution
		fileout.open("/Users/bpearce/octav/4/approxSol.txt");
		for (int i = 0; i < numsubintervals + 1; i++)
			fileout << sol[i] << endl;
		fileout.close();
        
        //If the true solution is present output it and the error in a few forms
		if (prob->true_solution_is_present())
		{
            //create a vector to save the true solution.
            vector<double>truSoln(numsubintervals+1);
            
			fileout.open("/Users/bpearce/octav/4/truesol.txt");

			double s = (dom[1] - dom[0]) / numsubintervals;
			vector<double> x(1);
			x[0] = dom[0];
			for (int i = 0; i < numsubintervals + 1; i++)
			{
                truSoln[i] = prob-> eval_true_solution(x);
                fileout << truSoln[i] << endl;
				x[0] += s;
			}
			fileout.close();
            
            //create and output an error vector for analysis
            vector <double> err(numsubintervals+1);
            fileout.open("/Users/bpearce/octav/4/approximationError.txt");
            for (int i = 0; i < numsubintervals+1; i++) {
                err[i] = truSoln[i] - sol[i];
                fileout << err[i] << endl;
            }
            fileout.close();
            
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
            
            double h = (dom[1] - dom[0])/numsubintervals;
            vector<double> deriv2(numsubintervals+1);
            deriv2 = scndDerv->Mult(sol);
            
            deriv2[0] = (sol[0] - 2*(sol[1]) + sol[2])/h;
            deriv2[numsubintervals] =(sol[numsubintervals-2]
                            - 2*(sol[numsubintervals-1])
                            + sol[numsubintervals])/h;
            
            // Output the 2nd Deriv to a output file
            fileout.open("/Users/bpearce/octav/4/scndDeriv.txt");
            for (int j=0; j < numsubintervals+1; j++) {
                fileout << deriv2[j] <<endl;
            }
            fileout.close();
		}
	}

	//----------------------------------------
	//end of section that finds the approx solution
	//----------------------------------------



	//-----------------------------------------
	//find error between true and approximate soln
	//(comment out if only finding approx soln
	// higlight then ctrl+k, ctrl+c to comment out
	// highlight then ctrl+k, ctrl+u to uncomment)
	//-----------------------------------------
	
 
  {
  
		{
			//create vectors to store hs and e(x_j)s and ln
			vector<double> h(10);
			vector<double> ln_h(10);
            vector<double> errors(2);
            vector<double> errL2(10);
			vector<double> ex_j(10);
			vector<double> ln_ex_j(10);
            vector<double> ln_L2_err(10);

			// for loop to run with diffrent size hs

			//counter to help update h & ex_j
			vector<double> doubles(10);
			//named doubles because each entry is a double of the previous
			//and this vector doubles the numsubintervals each iteration

            doubles[0] = 10;
            
			for (int i = 1; i < 10; i++)
			{
                doubles[i] = doubles[i-1]*2;
			}

			//for loop to run the appx over and over w/ diff h.
			for (int j = 0; j < 10; j++)
			{

				// create the 2 pt bvp approximation
				double numsubintervals = doubles[j];
				double * subintervals = new double[numsubintervals];
				for (int i = 0; i < numsubintervals; i++)
				{
					subintervals[i] = (dom[1] - dom[0]) / numsubintervals;
				}
				//create subintervals size vector
				h[j] = (dom[1] - dom[0]) / numsubintervals;

				// create the approximation for the new stepinterval size
				TwoPointBVPAppr *method = new TwoPointBVPAppr(numsubintervals,
					subintervals, prob);

				method->set_intial_guess_seed(seed);

				//solve and find the error vector and the max error
				//for the new step size h
                errors = method->find_errors(maxIters, Toler);
                errL2[j] = errors[0];
                ex_j[j] = errors[1];
			}

			//find the natural log of the subinterval lengths and errors
			//for plotting

			for (int i = 0; i < 10; i++)
			{
				ln_h[i] = log(h[i]);
				ln_ex_j[i] = log(ex_j[i]);
                ln_L2_err[i] = log(errL2[i]);
			}

			// output the subinterval lenght
			ofstream fileout;
			fileout.open("/Users/bpearce/octav/4/subintervallenght.txt");
			for (int i = 0; i < 10; i++)
			{
				fileout << h[i] << endl;
			}
			fileout.close();

            //output the error L2 norm for the cooresponding subinterval length
            fileout.open("/Users/bpearce/octav/4/L2error.txt");
            for (int i = 0; i < 10; i++)
            {
                fileout << errL2[i] <<  endl;
            }
            fileout.close();
            
            //output the max error for the cooresponding subinterval length
            fileout.open("/Users/bpearce/octav/4/max_error.txt");
            for (int i = 0; i < 10; i++)
            {
                fileout << ex_j[i] <<  endl;
            }
            fileout.close();
            
			// output the ln of subinterval lenght
			fileout.open("/Users/bpearce/octav/4/ln_subintervallenght.txt");
			for (int i = 0; i < 10; i++)
			{
				fileout << ln_h[i] << "\t" << ln_ex_j[i] << " " << endl;
			}
			fileout.close();
            
            // output the ln of the maximum error
            fileout.open("/Users/bpearce/octav/4/ln_max_error.txt");
            for (int i = 0; i < 10; i++)
            {
                fileout << ln_ex_j[i] << " " << endl;
            }
            fileout.close();
            
            fileout.open("/Users/bpearce/octav/4/ln_L2_error.txt");
            for (int i = 0; i < 10; i++)
            {
                fileout << ln_h[i] << "\t" << ln_ex_j[i] << " " << endl;
            }
            fileout.close();
			//delete vars associated with the problem
			//delete[]rbval;
			delete[]lbval;
			delete prob;
			delete[]dom;
		}
        
	}
	//------------------------------------------------------
	// end of section that finds l infinty norm
	//-------------------------------------------------------
  
	return 0;
}
