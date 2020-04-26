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
    std::mt19937 rng;
    rng.seed(std::random_device()());
    std::uniform_real_distribution<double> distribution(0, 1);
    
    std::uniform_int_distribution<int> distribution_uni(1,m_loci);
    double rndn=0;
    for (int i=1;i<=m_number;i++)
    {
        int mutation_loci=distribution_uni(rng);
        mutation_site(i,1)=mutation_location[mutation_loci-1];
        mutation_temp(i,1)=mutation_site(i,1);////////////////////
        rndn=distribution(rng);
        if(seq(mutation_site(i,1),1)==1)
        {
            if (rndn<sitefreq(4,2))
            {
                mutation_temp(i,2)=sitefreq(2,1);
            }
            else if (rndn>sitefreq(4,2) && rndn<sitefreq(4,3))
            {
                mutation_temp(i,2)=sitefreq(2,2);
            }
            else if (rndn>sitefreq(4,3))
            {
                mutation_temp(i,2)=sitefreq(2,3);
            }
        }
        if(seq(mutation_site(i,1),1)==2)
        {
            if (rndn<sitefreq(4,5))
            {
                mutation_temp(i,2)=sitefreq(2,4);
            }
            else if (rndn>sitefreq(4,5) && rndn<sitefreq(4,6))
            {
                mutation_temp(i,2)=sitefreq(2,5);
            }
            else if (rndn>sitefreq(4,6))
            {
                mutation_temp(i,2)=sitefreq(2,6);
            }
        }
        if(seq(mutation_site(i,1),1)==3)
        {
            if (rndn<sitefreq(4,8))
            {
                mutation_temp(i,2)=sitefreq(2,7);
            }
            else if (rndn>sitefreq(4,8) && rndn<sitefreq(4,9))
            {
                mutation_temp(i,2)=sitefreq(2,8);
            }
            else if (rndn>sitefreq(4,9))
            {
                mutation_temp(i,2)=sitefreq(2,9);
            }
        }
        if(seq(mutation_site(i,1),1)==4)
        {
            if (rndn<sitefreq(4,11))
            {
                mutation_temp(i,2)=sitefreq(2,10);
            }
            else if (rndn>sitefreq(4,11) && rndn<sitefreq(4,12))
            {
                mutation_temp(i,2)=sitefreq(2,11);
            }
            else if (rndn>sitefreq(4,12))
            {
                mutation_temp(i,2)=sitefreq(2,12);
            }
        }
    }
}
#endif /* mutation_hpp */
