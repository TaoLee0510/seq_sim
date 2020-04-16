//
//  selection.hpp
//  seq_simu
//
//  Created by Taolee on 3/30/20.
//  Copyright © 2020 Taolee. All rights reserved.
//

#ifndef selection_hpp
#define selection_hpp

#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <ctime>
#include <blitz/blitz.h>
#include <blitz/array.h>
using namespace std;
using namespace blitz;

void selection(int m_number,Array<int, 2> seq,Array<int,2> &mutation_temp,Array<int, 2> codon)
{
    Range all = Range::all();
    const gsl_rng_type *T;
    gsl_rng *r;
    gsl_rng_env_setup();
    T = gsl_rng_ranlxs0;
    gsl_rng_default_seed = ((unsigned long)(time(NULL)));
    r = gsl_rng_alloc(T);
    Array<int,2> codon_temp(1,m_number,FortranArray<2>());
    codon_temp=0;
    
    for (int i=1;i<=m_number;i++)
    {
        codon_temp(1,i)=mutation_temp(1,i)%3;
        
        Array<int,2> temp1(1,3,FortranArray<2>());
        temp1=0;
        Array<int,2> temp2(1,3,FortranArray<2>());
        temp2=0;
        if (codon_temp(1,i)==0)
        {
            temp1(1,1)=seq(1,mutation_temp(1,i)-2);
            temp1(1,2)=seq(1,mutation_temp(1,i)-1);
            temp1(1,3)=seq(1,mutation_temp(1,i));
            temp2(1,1)=seq(1,mutation_temp(1,i)-2);
            temp2(1,2)=seq(1,mutation_temp(1,i)-1);
            temp2(1,3)=mutation_temp(2,i);
            int code_old=0;
            int code_new=0;
            for (int j=1; j<=62; j++)
            {
                if (temp1(1,1)==codon(j,1)&&temp1(1,2)==codon(j,2)&&temp1(1,3)==codon(j,3))
                {
                    code_old=codon(j,4);
                }
                if (temp2(1,1)==codon(j,1)&&temp2(1,2)==codon(j,2)&&temp2(1,3)==codon(j,3))
                {
                    code_new=codon(j,4);
                }
            }
            if ((code_new==0 && code_old==0)||(code_new>0 && code_old==0)||(code_new==0 && code_old>0))
            {
                if (gsl_rng_uniform(r)>0.05)
                {
                    mutation_temp(all,i)=0;
                }
            }
            
        }
        else
        {
            if (gsl_rng_uniform(r)>0.05)
            {
                mutation_temp(all,i)=0;
            }
        }
    }
    gsl_rng_free(r);
}


#endif /* selection_hpp */
