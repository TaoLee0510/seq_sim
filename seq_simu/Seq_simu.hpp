//
//  Seq_simu.hpp
//  seq_simu
//
//  Created by Taolee on 4/14/20.
//  Copyright © 2020 Taolee. All rights reserved.
//

#ifndef Seq_simu_hpp
#define Seq_simu_hpp

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
#include <string>
#include <blitz/blitz.h>
#include <blitz/array.h>
#include "read_files.hpp"
#include "random_poisson.hpp"
#include "random_uniform.hpp"
#include "mutation.hpp"
#include "selection.hpp"
#include "save_file.hpp"
#include "save_file_2.hpp"
#include "save_file_3.hpp"
#include "NS_judge.hpp"
#include <math.h>
#include "confidence_calculation.hpp"
#include "pre_evolution.hpp"
#include "evolution.hpp"
#include "similarity_tg.hpp"
#include "random_mutation.hpp"
#include "extra_results.hpp"
using namespace std;
using namespace blitz;
void Seq_simu (int seq_length, double probability_of_mutation_sar,int DDAAYY,int DAY_pa,int interval,int times, int D, int random_mutation_number,int Divergency_sampling_times,string evolving_sequence, string reference_sequence, string codon_file, string substitution_file,int output_evo_seq_file)
{
    Range all = Range::all();
    //////////////////////////////////////////////////// Arrays initiation/////////////////////////////////////////////
    Array<int, 2> seq_orig_sar(1,seq_length,FortranArray<2>());//seq_orig matrix
    seq_orig_sar=0;
    Array<int, 2> seq_orig_tg(1,seq_length,FortranArray<2>());//seq_orig matrix
    seq_orig_tg=0;
    Array<int, 2> seq_sar(1,seq_length,FortranArray<2>());//seq matrix
    seq_sar=0;
    Array<int, 2> seq_sar1(1,seq_length,FortranArray<2>());//seq matrix
    seq_sar1=0;
    Array<int, 2> seq_sar2(1,seq_length,FortranArray<2>());//seq matrix
    seq_sar2=0;
    Array<int, 2> codon(62,4,FortranArray<2>());//condon
    codon=0;
    Array<double, 2> sitefreq(12,4,FortranArray<2>());//site alter freq
    sitefreq=0;
    Array<int,2> mutation_number_sar(1,DDAAYY,FortranArray<2>());//mutation numbers matrix
    mutation_number_sar=0;
    Array<int,2> M_number(1,DDAAYY,FortranArray<2>()); //mutation numbers matrix
    M_number=0;
    Array<int, 2> seq_evo_sar(D+1,seq_length,FortranArray<2>());//seq_evo
    seq_evo_sar=0;
    Array<int, 2> seq_evo_sar1(D+1,seq_length,FortranArray<2>());//seq_evo
    seq_evo_sar1=0;
    Array<int, 2> seq_evo_sar2(D+1,seq_length,FortranArray<2>());//seq_evo
    seq_evo_sar2=0;
    Array<int, 2> seq_rand_sar(D+1,seq_length,FortranArray<2>());//seq_rand
    seq_rand_sar=0;
    
    Array<double, 3> SA_out_temp(D,7,Divergency_sampling_times,FortranArray<3>());
    SA_out_temp=0;
    Array<double, 3> All_out_temp(D,7,Divergency_sampling_times,FortranArray<3>());
    All_out_temp=0;
    
    Array<double, 3> SA_out(D,11,times,FortranArray<3>());
    SA_out=0;
    Array<double, 3> All_out(D,11,times,FortranArray<3>());
    All_out=0;
    Array<double, 3> Similarity(D,4,times,FortranArray<3>());
    Similarity=0;
    Array<double, 2> Similarity_tg(1,4,FortranArray<2>());
    Similarity_tg=0;
    Array<double, 3> Similarity_divergency(D,4,times,FortranArray<3>());
    Similarity_divergency=0;
    
    
    //////////////////////////////////////////////////// read file /////////////////////////////////////////////
    read_file(seq_orig_sar, seq_orig_tg,codon, sitefreq,seq_length, evolving_sequence, reference_sequence, codon_file, substitution_file);
    similarity_tg(codon, seq_orig_sar, seq_orig_tg, Similarity_tg, seq_length,1);
    save_file_3(Similarity_tg);
    //////////////////////////////////////////////////// mutation location stastics /////////////////////////////////////////////
    int m_loci_sar=0;
    for (int i=1;i<=seq_length;i++)
    {
        if (seq_orig_sar(1,i)!=0)
        {
            m_loci_sar=m_loci_sar+1;
        }
    }
    int *mutation_location_sar=new int[m_loci_sar];
    int m_loci=0;
    for (int i=1;i<=seq_length;i++)
    {
        if (seq_orig_sar(1,i)!=0)
        {
            mutation_location_sar[m_loci]=i;
            m_loci=m_loci+1;
        }
    }
    m_loci=0;
    /////////////////////////////////////////////////// innitiation /////////////////////////////////////////////
    for (int duplic=1;duplic<=times;duplic++)
    {
        seq_evo_sar(1,all)=seq_orig_sar(1,all);
        seq_sar(1,all)=seq_orig_sar(1,all);
        seq_evo_sar1(1,all)=seq_orig_sar(1,all);
        seq_sar1(1,all)=seq_orig_sar(1,all);
        seq_sar2(1,all)=seq_orig_sar(1,all);
        //////////////////////////////////////////////////// pre evolution /////////////////////////////////////////////
        Array<int,2> mutation_number_sar1(1,DAY_pa,FortranArray<2>());//mutation numbers matrix
        mutation_number_sar1=0;
        mutation_number_sar1=random_poisson(DAY_pa, probability_of_mutation_sar);
        for (int day=1;day<=DAY_pa;day++)
        {
            pre_evolution(sitefreq, codon, seq_orig_sar, seq_sar2, mutation_number_sar1, m_loci_sar, mutation_location_sar, seq_length, day,duplic);
        }
        seq_evo_sar2(1,all)=seq_sar2(1,all);
        //////////////////////////////////////////////////// mutation number initiation /////////////////////////////////////////////
        mutation_number_sar=random_poisson(DDAAYY, probability_of_mutation_sar);
        //////////////////////////////////////////////////// evolution /////////////////////////////////////////////
        for (int day=1;day<=DDAAYY;day++)
        {
            evolution (sitefreq, codon, seq_orig_sar, seq_sar, seq_sar1, seq_sar2, seq_evo_sar, seq_evo_sar1, seq_evo_sar2, mutation_number_sar, m_loci_sar, mutation_location_sar, seq_length, interval, day, duplic);
        }
        if (output_evo_seq_file==1)
        {
            save_file(seq_evo_sar,seq_evo_sar1,seq_evo_sar2,duplic);
        }
        for (int rand_times=1;rand_times<=Divergency_sampling_times;rand_times++)
        {
            seq_rand_sar=0;
            seq_rand_sar(all,all)=seq_evo_sar(all,all);
            
            for (int j=2;j<=D+1;j++)
            {
                //////////////////////////////////////////////////// random mutation  /////////////////////////////////////////////
                random_mutation(sitefreq, codon, seq_rand_sar, m_loci_sar,mutation_location_sar, seq_length, j, random_mutation_number,duplic,rand_times);
                //////////////////////////////////////////////////// confidence  /////////////////////////////////////////////
                confidence_calculation(codon, seq_orig_sar, seq_evo_sar, seq_evo_sar1, seq_evo_sar2, seq_rand_sar, SA_out_temp, All_out_temp, Similarity_divergency, Similarity, j, rand_times, seq_length, duplic);
            }
        }
        extral_results(SA_out_temp, All_out_temp, SA_out, All_out, D, duplic);
        save_file_2(SA_out,All_out,Similarity,Similarity_divergency,D,duplic);
    }
    delete[] mutation_location_sar;
    mutation_location_sar = NULL;
}

#endif /* Seq_simu_hpp */
