//
//  Seq_simu.hpp
//  seq_simu
//
//  Created by Taolee on 3/29/20.
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
#include "accuracy.hpp"
#include "pre_evolution.hpp"
#include "evolution.hpp"
#include "forward_evolution.hpp"
#include "similarity_tg.hpp"
#include "extra_results.hpp"
#include <random>
using namespace std;
using namespace blitz;
void Seq_simu (int seq_length, double probability_of_mutation_sar, double mutation_rate_outgroup1, double mutation_rate_outgroup2, int DDAAYY,int DAY_pa,int interval,int times, int D, int random_mutation_number,int Divergency_sampling_times,string evolving_sequence, string reference_sequence, string codon_file, string substitution_file, string substitution_file_outgroup1, string substitution_file_outgroup2, int output_evo_seq_file,double KaKS,double Divergence, double n_selection)
{
    Range all = Range::all();
    //////////////////////////////////////////////////// Arrays initiation/////////////////////////////////////////////
    Array<int, 2> seq_orig_sar(seq_length,1,FortranArray<2>());//seq_orig matrix
    seq_orig_sar=0;
    Array<int, 2> seq_orig_tg(seq_length,1,FortranArray<2>());//seq_orig matrix
    seq_orig_tg=0;
    Array<int, 2> seq_sar(seq_length,1,FortranArray<2>());//seq matrix
    seq_sar=0;
    Array<int, 2> seq_sar1(seq_length,1,FortranArray<2>());//seq matrix
    seq_sar1=0;
    Array<int, 2> seq_sar2(seq_length,1,FortranArray<2>());//seq matrix
    seq_sar2=0;
    Array<int, 2> codon(4,62,FortranArray<2>());//condon
    codon=0;
    Array<double, 2> sitefreq(4,12,FortranArray<2>());//site alter freq
    sitefreq=0;
    Array<double, 2> sitefreq_out1(4,12,FortranArray<2>());//site alter freq
    sitefreq_out1=0;
    Array<double, 2> sitefreq_out2(4,12,FortranArray<2>());//site alter freq
    sitefreq_out2=0;
    Array<int,2> M_number(DDAAYY,1,FortranArray<2>()); //mutation numbers matrix
    M_number=0;
    Array<int, 2> seq_evo_sar(seq_length,D+1,FortranArray<2>());//seq_evo
    seq_evo_sar=0;
    Array<int, 2> seq_evo_sar1(seq_length,D+1,FortranArray<2>());//seq_evo
    seq_evo_sar1=0;
    Array<int, 2> seq_evo_sar2(seq_length,D+1,FortranArray<2>());//seq_evo
    seq_evo_sar2=0;
    Array<int, 2> seq_rand_sar(seq_length,D+1,FortranArray<2>());//seq_rand
    seq_rand_sar=0;
    Array<int, 2> seq_rand_sar1(seq_length,D+1,FortranArray<2>());//seq_rand
    seq_rand_sar1=0;
    Array<int, 2> seq_rand_sar0(seq_length,D+1,FortranArray<2>());//seq_rand
    seq_rand_sar0=0;
    Array<int, 2> seq_rand_sar01(seq_length,D+1,FortranArray<2>());//seq_rand
    seq_rand_sar01=0;
    Array<double, 3> SA_out_temp(9,D,Divergency_sampling_times,FortranArray<3>());
    SA_out_temp=0;
    Array<double, 3> All_out_temp(9,D,Divergency_sampling_times,FortranArray<3>());
    All_out_temp=0;
    Array<double, 3> SA_out(15,D,times,FortranArray<3>());
    SA_out=0;
    Array<double, 3> All_out(15,D,times,FortranArray<3>());
    All_out=0;
    Array<double, 3> Similarity(4,D,times,FortranArray<3>());
    Similarity=0;
    Array<double, 3> Similarity1(4,D,times,FortranArray<3>());
    Similarity1=0;
    Array<double, 2> Similarity_tg(1,4,FortranArray<2>());
    Similarity_tg=0;
    Array<double, 3> Similarity_divergency(4,D,times,FortranArray<3>());
    Similarity_divergency=0;
    Array<double, 3> Similarity_divergency1(4,D,times,FortranArray<3>());
    Similarity_divergency1=0;
    Array<double, 2> mutatioin_number(times,10,times,FortranArray<2>());
    mutatioin_number=0;
    //////////////////////////////////////////////////// read file /////////////////////////////////////////////
    read_file(seq_orig_sar, seq_orig_tg,codon, sitefreq,sitefreq_out1,sitefreq_out2, seq_length, evolving_sequence, reference_sequence, codon_file, substitution_file, substitution_file_outgroup1, substitution_file_outgroup2);
    similarity_tg(codon, seq_orig_sar, seq_orig_tg, Similarity_tg, seq_length,1);
    //////////////////////////////////////////////////// mutation location stastics /////////////////////////////////////////////
    int m_loci_sar=0;
    for (int i=1;i<=seq_length;i++)
    {
        if (seq_orig_sar(i,1)!=0)
        {
            m_loci_sar=m_loci_sar+1;
        }
    }
    int *mutation_location_sar=new int[m_loci_sar];
    int m_loci=0;
    for (int i=1;i<=seq_length;i++)
    {
        if (seq_orig_sar(i,1)!=0)
        {
            mutation_location_sar[m_loci]=i;
            m_loci=m_loci+1;
        }
    }
    m_loci=0;
    /////////////////////////////////////////////////// innitiation /////////////////////////////////////////////
    double selections=KaKS;
    double NS_site_number=seq_length*3.5/4.5;
    double S_site_number=seq_length/4.5;
    
    int forward_evo_time=0;
    if (n_selection>0)
    {
        forward_evo_time=ceil(((Divergence*seq_length)/(probability_of_mutation_sar*(1+(3.5*n_selection))/4.5))/2);
    }
    else
    {
        forward_evo_time=ceil((Divergence*m_loci_sar)/(probability_of_mutation_sar))/2;
    }
    cout << "forward_evo_time:"   << forward_evo_time << endl;
    int duplic=1;
    for (duplic=1;duplic<=times;duplic++)
    {
        seq_evo_sar(all,1)=seq_orig_sar(all,1);
        seq_sar(all,1)=seq_orig_sar(all,1);
        seq_evo_sar1(all,1)=seq_orig_sar(all,1);
        seq_sar1(all,1)=seq_orig_sar(all,1);
        seq_sar2(all,1)=seq_orig_sar(all,1);
        //////////////////////////////////////////////////// pre evolution /////////////////////////////////////////////
        int *mutation_number_sar0=new int[DAY_pa];
        mutation_number_sar0=random_poisson(DAY_pa, mutation_rate_outgroup2);
        int day=1;
        int mutation_real_s=0;
        int mutation_real_ns=0;
        for (day=1;day<=DAY_pa;day++)
        {
            pre_evolution(sitefreq_out2, codon, seq_orig_sar, seq_sar2, mutation_number_sar0, m_loci_sar, mutation_location_sar, seq_length, day,duplic,selections, mutation_real_s, mutation_real_ns);
        }
        seq_evo_sar2(all,1)=seq_sar2(all,1);
        //////////////////////////////////////////////////// mutation number initiation /////////////////////////////////////////////
        int *mutation_number_sar=new int[DDAAYY];
        mutation_number_sar=random_poisson(DDAAYY, probability_of_mutation_sar);
        int *mutation_number_sar1=new int[DDAAYY];
        mutation_number_sar1=random_poisson(DDAAYY, mutation_rate_outgroup1);
        int *mutation_number_sar2=new int[DDAAYY];
        mutation_number_sar2=random_poisson(DDAAYY, mutation_rate_outgroup2);
        //////////////////////////////////////////////////// evolution /////////////////////////////////////////////
        day=1;
        int total_mutation=0;
        int mutation_real=0;
        mutation_real_s=0;
        mutation_real_ns=0;
        for (day=1;day<=DDAAYY;day++)
        {
            evolution(sitefreq, sitefreq_out1, sitefreq_out2, codon, seq_orig_sar, seq_sar, seq_sar1, seq_sar2, seq_evo_sar, seq_evo_sar1, seq_evo_sar2, mutation_number_sar, mutation_number_sar1, mutation_number_sar2, m_loci_sar, mutation_location_sar, seq_length, interval, day, duplic, selections, mutation_real, mutation_real_s, mutation_real_ns);
            total_mutation=total_mutation+mutation_number_sar[day-1];
        }
        if (output_evo_seq_file==1)
        {
            save_file(seq_evo_sar,seq_evo_sar1,seq_evo_sar2,duplic);
        }
        int rand_times=1;
        int j=2;
        for (rand_times=1;rand_times<=Divergency_sampling_times;rand_times++)
        {
            cout <<"Runing times: "<<duplic<<"     Accuracy run times:    "<< rand_times << endl;
            seq_rand_sar=0;
            seq_rand_sar1=0;
            seq_rand_sar0=0;
            seq_rand_sar01=0;
            seq_rand_sar(all,all)=seq_evo_sar(all,all);
            seq_rand_sar1(all,all)=seq_evo_sar(all,all);
            seq_rand_sar0(all,all)=seq_evo_sar1(all,all);
            seq_rand_sar01(all,all)=seq_evo_sar2(all,all);
            for (j=2;j<=D+1;j++)
            {
                int *mutation_number_rand_sar=new int[forward_evo_time];
                mutation_number_rand_sar=random_poisson(forward_evo_time, probability_of_mutation_sar);
                int *mutation_number_rand_sar1=new int[forward_evo_time];
                mutation_number_rand_sar1=random_poisson(forward_evo_time, probability_of_mutation_sar);
                int *mutation_number_rand_sar0=new int[forward_evo_time];
                mutation_number_rand_sar0=random_poisson(forward_evo_time, mutation_rate_outgroup1);
                int *mutation_number_rand_sar01=new int[forward_evo_time];
                mutation_number_rand_sar01=random_poisson(forward_evo_time, mutation_rate_outgroup2);
                //////////////////////////////////////////////////// divergence  /////////////////////////////////////////////
                forward_evolution(sitefreq, sitefreq_out1, sitefreq_out2, codon, seq_rand_sar, seq_rand_sar1, seq_rand_sar0, seq_rand_sar01, mutation_number_rand_sar0, mutation_number_rand_sar, mutation_number_rand_sar1, mutation_number_rand_sar01, forward_evo_time, seq_length, j, m_loci_sar, selections,mutation_location_sar, n_selection);
                delete[] mutation_number_rand_sar;
                mutation_number_rand_sar = NULL;
                delete[] mutation_number_rand_sar1;
                mutation_number_rand_sar1 = NULL;
                delete[] mutation_number_rand_sar0;
                mutation_number_rand_sar0 = NULL;
                delete[] mutation_number_rand_sar01;
                mutation_number_rand_sar01 = NULL;
                //////////////////////////////////////////////////// accuracy  /////////////////////////////////////////////
                accuracy(codon, seq_orig_sar, seq_evo_sar, seq_evo_sar1, seq_evo_sar2, seq_rand_sar, seq_rand_sar1, seq_rand_sar0, seq_rand_sar01, SA_out_temp, All_out_temp, Similarity_divergency, Similarity_divergency1, Similarity, Similarity1, j, rand_times, seq_length, duplic, NS_site_number, S_site_number);
            }
            j=2;
        }
        extral_results(SA_out_temp, All_out_temp, SA_out, All_out, D, duplic);
        save_file_2(SA_out,All_out,Similarity,Similarity_divergency,Similarity_divergency1,Similarity1, D,duplic);
        mutatioin_number(duplic,1)=total_mutation;
        mutatioin_number(duplic,2)=mutation_real;
        mutatioin_number(duplic,3)=mutation_real/4.5;
        mutatioin_number(duplic,4)=(mutation_real/4.5)*0.175;
        mutatioin_number(duplic,5)=mutation_real;
        mutatioin_number(duplic,6)=mutation_real_s;
        mutatioin_number(duplic,7)=mutation_real_ns;
        mutatioin_number(duplic,8)=Similarity1(2,D,duplic)*seq_length;
        mutatioin_number(duplic,9)=Similarity1(3,D,duplic)*S_site_number;
        mutatioin_number(duplic,10)=Similarity1(4,D,duplic)*NS_site_number;
        delete[] mutation_number_sar0;
        mutation_number_sar0 = NULL;
        delete[] mutation_number_sar;
        mutation_number_sar = NULL;
        delete[] mutation_number_sar1;
        mutation_number_sar1 = NULL;
        delete[] mutation_number_sar2;
        mutation_number_sar2 = NULL;
    }
    save_file_3(Similarity_tg, mutatioin_number, times);
    delete[] mutation_location_sar;
    mutation_location_sar = NULL;
}
#endif /* Seq_simu_hpp */
