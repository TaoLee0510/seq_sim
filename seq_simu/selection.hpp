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
#include <random>
//#include "random_uniform.hpp"
using namespace std;
using namespace blitz;
void selection(int m_number,Array<int, 2> seq,Array<int,2> &mutation_temp,Array<int, 2> codon,double selections)
{
    Range all = Range::all();
    std::mt19937 rng;
    rng.seed(std::random_device()());
    std::uniform_real_distribution<double> distribution(0, 1);
    Array<int,2> temp1(3,1,FortranArray<2>());
    temp1=0;
    Array<int,2> temp2(3,1,FortranArray<2>());
    temp2=0;
    for (int i=1;i<=m_number;i++)
    {
        temp1=0;
        temp2=0;
        if (mutation_temp(i,1)%3==0)
        {
            temp1(1,1)=seq(mutation_temp(i,1)-2,1);
            temp1(2,1)=seq(mutation_temp(i,1)-1,1);
            temp1(3,1)=seq(mutation_temp(i,1),1);
            temp2(1,1)=seq(mutation_temp(i,1)-2,1);
            temp2(2,1)=seq(mutation_temp(i,1)-1,1);
            temp2(3,1)=mutation_temp(i,2);
            int code_old=1;
            int code_new=0;
            for (int j=1; j<=62; j++)
            {
                if (temp1(1,1)==codon(1,j)&&temp1(2,1)==codon(2,j)&&temp1(3,1)==codon(3,j))
                {
                    code_old=codon(4,j);
                }
                if (temp2(1,1)==codon(1,j)&&temp2(2,1)==codon(2,j)&&temp2(3,1)==codon(3,j))
                {
                    code_new=codon(4,j);
                }
            }
            if (code_new!=code_old)
            {
                if (distribution(rng)>selections)
                {
                    mutation_temp(i,all)=0;
                }
            }
        }
        else
        {
            if (distribution(rng)>selections)
            {
                mutation_temp(i,all)=0;
            }
        }
    }
}
#endif /* selection_hpp */
