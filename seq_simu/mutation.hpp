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
        mutation_site(i,1)=mutation_location[i-1];
        double *rndn=new double[1000];
        for (int i=0; i<1000; i++)
        {
            rndn[i]=gsl_rng_uniform(r1);
        }
        shuffle(rndn, rndn+1000,RNG);
        if(seq(mutation_site(i,1),1)==1)
        {
            if (rndn[0]<sitefreq(4,2))
            {
                mutation_temp(i,2)=sitefreq(2,1);
            }
            else if (rndn[0]>sitefreq(4,2) && rndn[0]<sitefreq(4,3))
            {
                mutation_temp(i,2)=sitefreq(2,2);
            }
            else if (rndn[0]>sitefreq(4,3))
            {
                mutation_temp(i,2)=sitefreq(2,3);
            }
        }
        if(seq(mutation_site(i,1),1)==2)
        {
            if (rndn[0]<sitefreq(4,5))
            {
                mutation_temp(i,2)=sitefreq(2,4);
            }
            else if (rndn[0]>sitefreq(4,5) && rndn[0]<sitefreq(4,6))
            {
                mutation_temp(i,2)=sitefreq(2,5);
            }
            else if (rndn[0]>sitefreq(4,6))
            {
                mutation_temp(i,2)=sitefreq(2,6);
            }
        }
        if(seq(mutation_site(i,1),1)==3)
        {
            if (rndn[0]<sitefreq(4,8))
            {
                mutation_temp(i,2)=sitefreq(2,7);
            }
            else if (rndn[0]>sitefreq(4,8) && rndn[0]<sitefreq(4,9))
            {
                mutation_temp(i,2)=sitefreq(2,8);
            }
            else if (rndn[0]>sitefreq(4,9))
            {
                mutation_temp(i,2)=sitefreq(2,9);
            }
        }
        if(seq(mutation_site(i,1),1)==4)
        {
            if (rndn[0]<sitefreq(4,11))
            {
                mutation_temp(i,2)=sitefreq(2,10);
            }
            else if (rndn[0]>sitefreq(4,11) && rndn[0]<sitefreq(4,12))
            {
                mutation_temp(i,2)=sitefreq(2,11);
            }
            else if (rndn[0]>sitefreq(4,12))
            {
                mutation_temp(i,2)=sitefreq(2,12);
            }
        }
        delete[] rndn;
        rndn= NULL;
    }
    mutation_temp(all,1)=mutation_site(all,1);
    gsl_rng_free(r1);
}
#endif /* mutation_hpp */
