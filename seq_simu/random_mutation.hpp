//
//  random_mutation.hpp
//  seq_simu
//
//  Created by Taolee on 3/29/20.
//  Copyright © 2020 Taolee. All rights reserved.
//

#ifndef random_mutation_hpp
#define random_mutation_hpp

#include <iostream>
#include <time.h>
#include <memory>
#include <stdio.h>
#include <cmath>
#include <algorithm>
#include <functional>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <getopt.h>
#include <blitz/blitz.h>
#include <blitz/array.h>
#include "random_poisson.hpp"
#include "mutation.hpp"
#include "selection.hpp"
#include <math.h>

using namespace std;
using namespace blitz;
void random_mutation (Array<double, 2> sitefreq, Array<int, 2> codon,Array<int, 2> &seq_rand_sar,int m_loci_sar,int *mutation_location_sar,int seq_length,int j,int random_mutation_number,int duplic,int rand_times,double selections)
{
    Range all = Range::all();
    Array<int,2> seq_sar(seq_length,1,FortranArray<2>());
    seq_sar=0;
    seq_sar(all,1)=seq_rand_sar(all,j);
    Array<int,2> mutation_site(random_mutation_number,1,FortranArray<2>());
    mutation_site=0;
    Array<int,2> mutation_temp(random_mutation_number,2,FortranArray<2>());
    mutation_temp=0;
    mutation(random_mutation_number, seq_sar, mutation_site, mutation_temp, sitefreq,seq_length,mutation_location_sar,m_loci_sar);//mutation
    selection(random_mutation_number, seq_sar, mutation_temp,codon,selections);//selection
    for (int i=1; i<=random_mutation_number;i++)
    {
        if (mutation_temp(i,1)>0)
        {
            seq_sar(mutation_temp(i,1),1)=mutation_temp(i,2);
        }
    }
    seq_rand_sar(all,j)= seq_sar(all,1);
}
#endif /* random_mutation_hpp */
