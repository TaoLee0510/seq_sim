//
//  forward_evolution.hpp
//  seq_simu
//
//  Created by Taolee on 4/24/20.
//  Copyright © 2020 Taolee. All rights reserved.
//

#ifndef forward_evolution_hpp
#define forward_evolution_hpp

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

void forward_evolution (Array<double, 2> sitefreq, Array<double, 2> sitefreq_out1,Array<double, 2> sitefreq_out2,Array<int, 2> codon, Array<int, 2> &seq_rand_sar, Array<int, 2> &seq_rand_sar1, Array<int, 2> &seq_rand_sar0, Array<int, 2> &seq_rand_sar01, int *mutation_number_rand_sar0, int *mutation_number_rand_sar, int *mutation_number_rand_sar1, int *mutation_number_rand_sar01, int forward_evo_time, int seq_length, int j, int m_loci_sar, double selections,int *mutation_location_sar)
{
    Range all = Range::all();
    Array<int,2> seq_sar(seq_length,1,FortranArray<2>());
    seq_sar=0;
    seq_sar(all,1)=seq_rand_sar(all,j);
    
    Array<int,2> seq_sar1(seq_length,1,FortranArray<2>());
    seq_sar1=0;
    seq_sar1(all,1)=seq_rand_sar1(all,j);
    
    Array<int,2> seq_sar0(seq_length,1,FortranArray<2>());
    seq_sar0=0;
    seq_sar0(all,1)=seq_rand_sar0(all,j);
    
    Array<int,2> seq_sar01(seq_length,1,FortranArray<2>());
    seq_sar01=0;
    seq_sar01(all,1)=seq_rand_sar01(all,j);
    
    for (int t3=1;t3<=forward_evo_time;t3++)
    {
        int temp1=0;
        int temp2=0;
        
        int m_number_sar=mutation_number_rand_sar[t3-1];
        if (m_number_sar>0)
        {
            Array<int,2> mutation_site(m_number_sar,1,FortranArray<2>());
            mutation_site=0;
            Array<int,2> mutation_temp(m_number_sar,2,FortranArray<2>());
            mutation_temp=0;
            mutation(m_number_sar, seq_sar, mutation_site, mutation_temp, sitefreq,seq_length,mutation_location_sar,m_loci_sar);
            selection(m_number_sar, seq_sar, mutation_temp,codon,selections,temp1,temp2);
            for (int i=1; i<=m_number_sar;i++)
            {
                if (mutation_temp(i,1)>0)
                {
                    seq_sar(mutation_temp(i,1),1)=mutation_temp(i,2);
                }
            }
        }
        int m_number_sar1=mutation_number_rand_sar1[t3-1];
        if (m_number_sar1>0)
        {
            Array<int,2> mutation_site(m_number_sar1,1,FortranArray<2>());
            mutation_site=0;
            Array<int,2> mutation_temp(m_number_sar1,2,FortranArray<2>());
            mutation_temp=0;
            mutation(m_number_sar1, seq_sar1, mutation_site, mutation_temp, sitefreq,seq_length,mutation_location_sar,m_loci_sar);
            selection(m_number_sar1, seq_sar1, mutation_temp,codon,selections, temp1, temp2);
            for (int i=1; i<=m_number_sar1;i++)
            {
                if (mutation_temp(i,1)>0)
                {
                    seq_sar1(mutation_temp(i,1),1)=mutation_temp(i,2);
                }
            }
        }
        int m_number_sar0=mutation_number_rand_sar0[t3-1];
        if (m_number_sar0>0)
        {
            Array<int,2> mutation_site(m_number_sar0,1,FortranArray<2>());
            mutation_site=0;
            Array<int,2> mutation_temp(m_number_sar0,2,FortranArray<2>());
            mutation_temp=0;
            mutation(m_number_sar0, seq_sar0, mutation_site, mutation_temp, sitefreq_out1,seq_length,mutation_location_sar,m_loci_sar);
            selection(m_number_sar0, seq_sar0, mutation_temp,codon,selections, temp1, temp2);
            for (int i=1; i<=m_number_sar0;i++)
            {
                if (mutation_temp(i,1)>0)
                {
                    seq_sar0(mutation_temp(i,1),1)=mutation_temp(i,2);
                }
            }
        }
        int m_number_sar01=mutation_number_rand_sar01[t3-1];
        if (m_number_sar01>0)
        {
            Array<int,2> mutation_site(m_number_sar01,1,FortranArray<2>());
            mutation_site=0;
            Array<int,2> mutation_temp(m_number_sar01,2,FortranArray<2>());
            mutation_temp=0;
            mutation(m_number_sar01, seq_sar01, mutation_site, mutation_temp, sitefreq_out2,seq_length,mutation_location_sar,m_loci_sar);
            selection(m_number_sar01, seq_sar01, mutation_temp,codon,selections, temp1, temp2);
            for (int i=1; i<=m_number_sar01;i++)
            {
                if (mutation_temp(i,1)>0)
                {
                    seq_sar01(mutation_temp(i,1),1)=mutation_temp(i,2);
                }
            }
        }
    }
    seq_rand_sar(all,j)= seq_sar(all,1);
    seq_rand_sar1(all,j)= seq_sar1(all,1);
    seq_rand_sar0(all,j)= seq_sar0(all,1);
    seq_rand_sar01(all,j)= seq_sar01(all,1);
}

#endif /* forward_evolution_hpp */
