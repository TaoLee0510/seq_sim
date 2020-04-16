//
//  evolution.hpp
//  seq_simu
//
//  Created by Taolee on 3/29/20.
//  Copyright © 2020 Taolee. All rights reserved.
//

#ifndef evolution_hpp
#define evolution_hpp

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
void evolution (Array<double, 2> sitefreq, Array<int, 2> codon,Array<int, 2> seq_orig_sar,Array<int, 2> &seq_sar,Array<int, 2> &seq_sar1,Array<int, 2> &seq_sar2,Array<int, 2> &seq_evo_sar,Array<int, 2> &seq_evo_sar1,Array<int, 2> &seq_evo_sar2,Array<int,2> mutation_number_sar,int m_loci_sar,  int *mutation_location_sar,int seq_length,int interval,int day,int duplic,double KaKS)
{
    Range all = Range::all();
    int m_number_sar=mutation_number_sar(day,1);
    if (m_number_sar>0)
    {
        Array<int,2> mutation_site(m_number_sar,1,FortranArray<2>());
        mutation_site=0;
        Array<int,2> mutation_temp(m_number_sar,2,FortranArray<2>());
        mutation_temp=0;
        mutation(m_number_sar, seq_sar, mutation_site, mutation_temp, sitefreq,seq_length,mutation_location_sar,m_loci_sar);//mutation
        selection(m_number_sar, seq_sar, mutation_temp,codon,KaKS);//selection
        for (int i=1; i<=m_number_sar;i++)
        {
            if (mutation_temp(i,1)>0)
            {
                seq_sar(mutation_temp(i,1),1)=mutation_temp(i,2);
            }
        }
    }
    if (m_number_sar>0)
    {
        Array<int,2> mutation_site(m_number_sar,1,FortranArray<2>());
        mutation_site=0;
        Array<int,2> mutation_temp(m_number_sar,2,FortranArray<2>());
        mutation_temp=0;
        mutation(m_number_sar, seq_sar1, mutation_site, mutation_temp, sitefreq,seq_length,mutation_location_sar,m_loci_sar);//mutation
        selection(m_number_sar, seq_sar1, mutation_temp,codon,KaKS);//selection
        for (int i=1; i<=m_number_sar;i++)
        {
            if (mutation_temp(i,1)>0)
            {
                seq_sar1(mutation_temp(i,1),1)=mutation_temp(i,2);
            }
        }
    }
    
    if (m_number_sar>0)
    {
        Array<int,2> mutation_site(m_number_sar,1,FortranArray<2>());
        mutation_site=0;
        Array<int,2> mutation_temp(m_number_sar,2,FortranArray<2>());
        mutation_temp=0;
        mutation(m_number_sar, seq_sar2, mutation_site, mutation_temp, sitefreq,seq_length,mutation_location_sar,m_loci_sar);//mutation
        selection(m_number_sar, seq_sar2, mutation_temp,codon,KaKS);//selection
        for (int i=1; i<=m_number_sar;i++)
        {
            if (mutation_temp(i,1)>0)
            {
                seq_sar2(mutation_temp(i,1),1)=mutation_temp(i,2);
            }
        }
    }
    
    
    if (day%interval==0)
    {
        int D=day/interval;
        seq_evo_sar(all,D+1)=seq_sar(all,1);
        seq_evo_sar1(all,D+1)=seq_sar1(all,1);
        seq_evo_sar2(all,D+1)=seq_sar2(all,1);
        cout <<"Runing times: "<<duplic<<"     "<<"Day: "<< day << endl;
    }
}
#endif /* evolution_hpp */