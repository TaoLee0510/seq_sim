//
//  mutation.hpp
//  seq_simu
//
//  Created by Taolee on 3/30/20.
//  Copyright © 2020 Taolee. All rights reserved.
//

#ifndef mutation_hpp
#define mutation_hpp

#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <ctime>
#include <blitz/blitz.h>
#include <blitz/array.h>
#include <random>
using namespace blitz;
void mutation(int m_number,Array<int, 2> seq,Array<int,2> &mutation_site,Array<int,2> &mutation_temp, Array<double, 2> sitefreq,int seq_length,int *mutation_location, int m_loci)
{
    Range all = Range::all();
    
    std::random_device r;
    std::seed_seq seed{r(), r(), r(), r(), r(), r(), r(), r()};
    std::mt19937 RNG(seed);
    
    const gsl_rng_type *T1;
    gsl_rng *r1;
    gsl_rng_env_setup();
    T1 = gsl_rng_ranlxs0;
    gsl_rng_default_seed = ((unsigned long)(time(NULL)));
    r1 = gsl_rng_alloc(T1);
    shuffle(mutation_location, mutation_location+m_loci,RNG);
    for (int i=1;i<=m_number;i++)
    {
        mutation_site(1,i)=mutation_location[i-1];

        double *rndn=new double[1000];;
        for (int i=0; i<1000; i++)
        {
            rndn[i]=gsl_rng_uniform(r1);
        }
        shuffle(rndn, rndn+1000,RNG);
        if(seq(1,mutation_site(1,i))==1)
        {
            if (rndn[0]<sitefreq(2,4))
            {
                mutation_temp(2,i)=sitefreq(1,2);
            }
            else if (rndn[0]>sitefreq(2,4) && rndn[0]<sitefreq(3,4))
            {
                mutation_temp(2,i)=sitefreq(2,2);
            }
            else if (rndn[0]>sitefreq(3,4))
            {
                mutation_temp(2,i)=sitefreq(3,2);
            }
        }
        if(seq(1,mutation_site(1,i))==2)
        {
            if (rndn[0]<sitefreq(5,4))
            {
                mutation_temp(2,i)=sitefreq(4,2);
            }
            else if (rndn[0]>sitefreq(5,4) && rndn[0]<sitefreq(6,4))
            {
                mutation_temp(2,i)=sitefreq(5,2);
            }
            else if (rndn[0]>sitefreq(6,4))
            {
                mutation_temp(2,i)=sitefreq(6,2);
            }
        }
        if(seq(1,mutation_site(1,i))==3)
        {
            if (rndn[0]<sitefreq(8,4))
            {
                mutation_temp(2,i)=sitefreq(7,2);
            }
            else if (rndn[0]>sitefreq(8,4) && rndn[0]<sitefreq(9,4))
            {
                mutation_temp(2,i)=sitefreq(8,2);
            }
            else if (rndn[0]>sitefreq(9,4))
            {
                mutation_temp(2,i)=sitefreq(9,2);
            }
        }
        if(seq(1,mutation_site(1,i))==4)
        {
            if (rndn[0]<sitefreq(11,4))
            {
                mutation_temp(2,i)=sitefreq(10,2);
            }
            else if (rndn[0]>sitefreq(11,4) && rndn[0]<sitefreq(12,4))
            {
                mutation_temp(2,i)=sitefreq(11,2);
            }
            else if (rndn[0]>sitefreq(12,4))
            {
                mutation_temp(2,i)=sitefreq(12,2);
            }
        }
        delete[] rndn;
        rndn= NULL;
    }
    mutation_temp(1,all)=mutation_site(1,all);
    gsl_rng_free(r1);
}
#endif /* mutation_hpp */
