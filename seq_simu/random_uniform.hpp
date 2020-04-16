//
//  random_uniform.hpp
//  seq_simu
//
//  Created by Taolee on 3/30/20.
//  Copyright © 2020 Taolee. All rights reserved.
//

#ifndef random_uniform_hpp
#define random_uniform_hpp

#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <ctime>
#include <blitz/blitz.h>
#include <blitz/array.h>
using namespace blitz;
Array<float,2> random_uniform (int N0)
{
    const gsl_rng_type *T1;
    gsl_rng *r1;
    gsl_rng_env_setup();
    T1 = gsl_rng_ranlxs0;
    gsl_rng_default_seed = ((unsigned long)(time(NULL)));
    r1 = gsl_rng_alloc(T1);
    Array<float,2> A(N0,1,FortranArray<2>());
    for(int x=1;x<=N0;x++)
    {
        A(x,1) = gsl_rng_uniform(r1);
    }
    return A;
    gsl_rng_free(r1);
}

#endif /* random_uniform_hpp */
