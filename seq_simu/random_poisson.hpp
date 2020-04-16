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
#include <blitz/blitz.h>
#include <blitz/array.h>
using namespace blitz;
Array<int,2> random_poisson(int N0, double probability_of_mutation)
{
    const gsl_rng_type *T;
    gsl_rng *r;
    gsl_rng_env_setup();
    T = gsl_rng_ranlxs0;
    gsl_rng_default_seed = ((unsigned long)(time(NULL)));
    r = gsl_rng_alloc(T);
    Array<int,2> A(N0,1,FortranArray<2>());
    for(int x=1;x<=N0;x++)
    {
        A(x,1) = gsl_ran_poisson(r, probability_of_mutation);
    }
    return A;
    gsl_rng_free(r);
}

#endif /* random_poisson_hpp */
