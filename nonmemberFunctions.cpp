//
//  nonmemberFunctions.cpp
//  Torch
//
//  Created by Brent Morrill Pearce on 7/19/18.
//  Copyright © 2018 Brent Morrill Pearce. All rights reserved.
//

#include "nonmemberFucntions.h"

// Function to find the second norm.

double findL2Norm( vector<double> const x)
{
    int num_entries = x.size();
    double sum_o_squares = 0;
    for (int i = 0; i < num_entries ; i++)
    {
        sum_o_squares = sum_o_squares + x[i] * x[i];
    }
    return sqrt(sum_o_squares);
}
