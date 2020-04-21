//
//  random_poisson.hpp
//  seq_simu
//
//  Created by Taolee on 3/30/20.
//  Copyright © 2020 Taolee. All rights reserved.
//

#ifndef random_poisson_hpp
#define random_poisson_hpp

#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <ctime>
#include <random>
using namespace std;
int * random_poisson(int N0, double probability_of_mutation)
{
    std::mt19937 rng;
    rng.seed(std::random_device()());
    std::poisson_distribution<int> distribution(probability_of_mutation);
    int *A =new int[N0];
    for(int x=0;x<N0;x++)
    {
        A[x]=distribution(rng);
    }
    return A;
    delete[] A;
    A = NULL;
}

#endif /* random_poisson_hpp */
