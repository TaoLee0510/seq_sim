//
//  pre_evolution.hpp
//  seq_simu
//
//  Created by Taolee on 4/10/20.
//  Copyright © 2020 Taolee. All rights reserved.
//

#ifndef pre_evolution_hpp
#define pre_evolution_hpp

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
#include "random_uniform.hpp"
#include "mutation.hpp"
#include "selection.hpp"
#include <math.h>
using namespace std;
using namespace blitz;
void pre_evolution (Array<double, 2> sitefreq, Array<int, 2> codon,Array<int, 2> seq_orig_sar,Array<int, 2> &seq_sar2,Array<int,2> mutation_number_sar1,int m_loci_sar,int *mutation_location_sar,int seq_length, int day,int duplic)
{
    int m_number_sar=mutation_number_sar1(1,day);
    if (m_number_sar>0)
    {
        Array<int,2> mutation_site(1,m_number_sar,FortranArray<2>());
        mutation_site=0;
        Array<int,2> mutation_temp(2,m_number_sar,FortranArray<2>());
        mutation_temp=0;
        mutation(m_number_sar, seq_sar2, mutation_site, mutation_temp, sitefreq,seq_length,mutation_location_sar,m_loci_sar);//mutation
        selection(m_number_sar, seq_sar2, mutation_temp,codon);//selection
        for (int i=1; i<=m_number_sar;i++)
        {
            if (mutation_temp(1,i)>0)
            {
                seq_sar2(1,mutation_temp(1,i))=mutation_temp(2,i);
            }
        }
    }
    cout <<"pre evolution: "<<duplic<<"       Day: "<< day << endl;
}

#endif /* pre_evolution_hpp */
