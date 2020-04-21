//
//  confidence_calculation.hpp
//  seq_simu
//
//  Created by Taolee on 3/29/20.
//  Copyright © 2020 Taolee. All rights reserved.
//

#ifndef confidence_calculation_hpp
#define confidence_calculation_hpp

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
#include "NS_judge.hpp"
#include <math.h>

using namespace std;
using namespace blitz;
void confidence_calculation (Array<int, 2> codon, Array<int, 2> seq_orig_sar, Array<int, 2> seq_evo_sar,Array<int, 2> seq_evo_sar1,Array<int, 2> seq_evo_sar2, Array<int, 2> seq_rand_sar,Array<double, 3> &SA_out_temp, Array<double, 3> &All_out_temp ,Array<double, 3> &Similarity_divergency,Array<double, 3> &Similarity,int j, int rand_times, int seq_length,int duplic, double NS_site_number, double S_site_number)
{
    if (rand_times==1)
    {
        Range all = Range::all();
        int N_right_All_out=0;
        int N_err_All_out=0;
        int N_right_All_out_S=0;
        int N_err_All_out_S=0;
        int N_uncertainty_all=0;
        int N_right_SA_out=0;
        int N_err_SA_out=0;
        int N_right_SA_out_S=0;
        int N_err_SA_out_S=0;
        int N_uncertainty=0;
        int total_mutations=0;
        int NS_similarity=0;
        int S_similarity=0;
        int NS_similarity1=0;
        int S_similarity1=0;
        for (int i=1;i<=seq_length;i++)
        {
            if (seq_rand_sar(i,j)!=seq_evo_sar(i,j))
            {
                if ((seq_evo_sar1(i,j)==seq_evo_sar2(i,j) && seq_evo_sar1(i,j)==seq_evo_sar(i,j)) || (seq_evo_sar1(i,j)==seq_evo_sar(i,j) && seq_evo_sar2(i,j)!=seq_rand_sar(i,j))|| (seq_evo_sar2(i,j)==seq_evo_sar(i,j) && seq_evo_sar1(i,j)!=seq_rand_sar(i,j)))
                {
                    N_right_All_out=N_right_All_out+1;
                }
                else if ((seq_evo_sar1(i,j)==seq_evo_sar2(i,j) && seq_evo_sar1(i,j)==seq_rand_sar(i,j)) || (seq_evo_sar1(i,j)==seq_rand_sar(i,j) && seq_evo_sar2(i,j)!=seq_evo_sar(i,j)) || (seq_evo_sar2(i,j)==seq_rand_sar(i,j) && seq_evo_sar1(i,j)!=seq_evo_sar(i,j)))
                {
                    N_err_All_out=N_err_All_out+1;
                }
                else if ((seq_evo_sar2(i,j)!=seq_rand_sar(i,j) && seq_evo_sar2(i,j)!=seq_evo_sar(i,j)) && (seq_evo_sar1(i,j)!=seq_rand_sar(i,j) && seq_evo_sar1(i,j)!=seq_evo_sar(i,j)))
                {
                    N_uncertainty_all=N_uncertainty_all+1;
                }
                if (seq_evo_sar(i,j)==seq_evo_sar1(i,j))
                {
                    N_right_SA_out=N_right_SA_out+1;
                }
                else if (seq_evo_sar1(i,j)==seq_rand_sar(i,j))
                {
                    N_err_SA_out=N_err_SA_out+1;
                }
                else
                {
                    N_uncertainty=N_uncertainty+1;
                }
                total_mutations=total_mutations+1;
                if (i%3==0)
                {
                    Array<int, 2> seqsar(seq_length,1,FortranArray<2>());//seq_rand
                    seqsar=0;
                    seqsar(all,1)=seq_evo_sar(all,j);
                    int s=NS_judge(i,seqsar,seq_rand_sar,codon,j);
                    if(s==1)
                    {
                        if ((seq_evo_sar1(i,j)==seq_evo_sar2(i,j) && seq_evo_sar1(i,j)==seq_evo_sar(i,j)) || (seq_evo_sar1(i,j)==seq_evo_sar(i,j) && seq_evo_sar2(i,j)!=seq_rand_sar(i,j))|| (seq_evo_sar2(i,j)==seq_evo_sar(i,j) && seq_evo_sar1(i,j)!=seq_rand_sar(i,j)))
                        {
                            N_right_All_out_S=N_right_All_out_S+1;
                        }
                        else if ((seq_evo_sar1(i,j)==seq_evo_sar2(i,j) && seq_evo_sar1(i,j)==seq_rand_sar(i,j)) || (seq_evo_sar1(i,j)==seq_rand_sar(i,j) && seq_evo_sar2(i,j)!=seq_evo_sar(i,j)) || (seq_evo_sar2(i,j)==seq_rand_sar(i,j) && seq_evo_sar1(i,j)!=seq_evo_sar(i,j)))
                        {
                            N_err_All_out_S=N_err_All_out_S+1;
                        }
                        if (seq_evo_sar(i,j)==seq_evo_sar1(i,j))
                        {
                            N_right_SA_out_S=N_right_SA_out_S+1;
                        }
                        else if (seq_evo_sar1(i,j)==seq_rand_sar(i,j))
                        {
                            N_err_SA_out_S=N_err_SA_out_S+1;
                        }
                    }
                }
            }
            if (seq_evo_sar1(i,j)!=seq_evo_sar(i,j))
            {
                if (i%3==0)
                {
                    Array<int, 2> seqsar(seq_length,1,FortranArray<2>());//seq_rand
                    seqsar=0;
                    seqsar(all,1)=seq_evo_sar1(all,j);
                    int s=NS_judge(i,seqsar,seq_evo_sar,codon,j);
                    if (s==1)
                    {
                        S_similarity=S_similarity+1;
                    }
                    else
                    {
                        NS_similarity=NS_similarity+1;
                    }
                }
                else
                {
                    NS_similarity=NS_similarity+1;
                }
            }
            if (seq_rand_sar(i,j)!=seq_evo_sar(i,j))
            {
                if (i%3==0)
                {
                    Array<int, 2> seqsar(seq_length,1,FortranArray<2>());//seq_rand
                    seqsar=0;
                    seqsar(all,1)=seq_rand_sar(all,j);
                    int s=NS_judge(i,seqsar,seq_evo_sar,codon,j);
                    if (s==1)
                    {
                        S_similarity1=S_similarity1+1;
                    }
                    else
                    {
                        NS_similarity1=NS_similarity1+1;
                    }
                }
                else
                {
                    NS_similarity1=NS_similarity1+1;
                }
            }
        }
        SA_out_temp(1,j-1,rand_times)=N_right_SA_out;
        SA_out_temp(2,j-1,rand_times)=N_err_SA_out;
        SA_out_temp(3,j-1,rand_times)=N_right_SA_out_S;
        SA_out_temp(4,j-1,rand_times)=N_err_SA_out_S;
        SA_out_temp(5,j-1,rand_times)=N_right_SA_out-N_right_SA_out_S;
        SA_out_temp(6,j-1,rand_times)=N_err_SA_out-N_err_SA_out_S;
        SA_out_temp(7,j-1,rand_times)=N_uncertainty;
        All_out_temp(1,j-1,rand_times)=N_right_All_out;
        All_out_temp(2,j-1,rand_times)=N_err_All_out;
        All_out_temp(3,j-1,rand_times)=N_right_All_out_S;
        All_out_temp(4,j-1,rand_times)=N_err_All_out_S;
        All_out_temp(5,j-1,rand_times)=N_right_All_out-N_right_All_out_S;
        All_out_temp(6,j-1,rand_times)=N_err_All_out-N_err_All_out_S;
        All_out_temp(7,j-1,rand_times)=N_uncertainty_all;
        Similarity(1,j-1,duplic)=j;
        double sum5=S_similarity+NS_similarity;
        Similarity(2,j-1,duplic)=sum5/seq_length;
        Similarity(3,j-1,duplic)=S_similarity/S_site_number;
        Similarity(4,j-1,duplic)=NS_similarity/NS_site_number;
        Similarity_divergency(1,j-1,duplic)=j;
        double sum6=S_similarity1+NS_similarity1;
        Similarity_divergency(2,j-1,duplic)=sum6/seq_length;
        Similarity_divergency(3,j-1,duplic)=S_similarity1/S_site_number;
        Similarity_divergency(4,j-1,duplic)=NS_similarity1/NS_site_number;
    }
    else
    {
        Range all = Range::all();
        int N_right_All_out=0;
        int N_err_All_out=0;
        int N_right_All_out_S=0;
        int N_err_All_out_S=0;
        int N_uncertainty_all=0;
        int N_right_SA_out=0;
        int N_err_SA_out=0;
        int N_right_SA_out_S=0;
        int N_err_SA_out_S=0;
        int N_uncertainty=0;
        int total_mutations=0;
        for (int i=1;i<=seq_length;i++)
        {
            if (seq_rand_sar(i,j)!=seq_evo_sar(i,j))
            {
                if ((seq_evo_sar1(i,j)==seq_evo_sar2(i,j) && seq_evo_sar1(i,j)==seq_evo_sar(i,j)) || (seq_evo_sar1(i,j)==seq_evo_sar(i,j) && seq_evo_sar2(i,j)!=seq_rand_sar(i,j))|| (seq_evo_sar2(i,j)==seq_evo_sar(i,j) && seq_evo_sar1(i,j)!=seq_rand_sar(i,j)))
                {
                    N_right_All_out=N_right_All_out+1;
                }
                else if ((seq_evo_sar1(i,j)==seq_evo_sar2(i,j) && seq_evo_sar1(i,j)==seq_rand_sar(i,j)) || (seq_evo_sar1(i,j)==seq_rand_sar(i,j) && seq_evo_sar2(i,j)!=seq_evo_sar(i,j)) || (seq_evo_sar2(i,j)==seq_rand_sar(i,j) && seq_evo_sar1(i,j)!=seq_evo_sar(i,j)))
                {
                    N_err_All_out=N_err_All_out+1;
                }
                else if ((seq_evo_sar2(i,j)!=seq_rand_sar(i,j) && seq_evo_sar2(i,j)!=seq_evo_sar(i,j)) && (seq_evo_sar1(i,j)!=seq_rand_sar(i,j) && seq_evo_sar1(i,j)!=seq_evo_sar(i,j)))
                {
                    N_uncertainty_all=N_uncertainty_all+1;
                }
                if (seq_evo_sar(i,j)==seq_evo_sar1(i,j))
                {
                    N_right_SA_out=N_right_SA_out+1;
                }
                else if (seq_evo_sar1(i,j)==seq_rand_sar(i,j))
                {
                    N_err_SA_out=N_err_SA_out+1;
                }
                else
                {
                    N_uncertainty=N_uncertainty+1;
                }
                total_mutations=total_mutations+1;
                if (i%3==0)
                {
                   Array<int, 2> seqsar(seq_length,1,FortranArray<2>());//seq_rand
                    seqsar=0;
                    seqsar(all,1)=seq_evo_sar(all,j);
                    int s=NS_judge(i,seqsar,seq_rand_sar,codon,j);
                    if(s==1)
                    {
                        if ((seq_evo_sar1(i,j)==seq_evo_sar2(i,j) && seq_evo_sar1(i,j)==seq_evo_sar(i,j)) || (seq_evo_sar1(i,j)==seq_evo_sar(i,j) && seq_evo_sar2(i,j)!=seq_rand_sar(i,j))|| (seq_evo_sar2(i,j)==seq_evo_sar(i,j) && seq_evo_sar1(i,j)!=seq_rand_sar(i,j)))
                        {
                            N_right_All_out_S=N_right_All_out_S+1;
                        }
                        else if ((seq_evo_sar1(i,j)==seq_evo_sar2(i,j) && seq_evo_sar1(i,j)==seq_rand_sar(i,j)) || (seq_evo_sar1(i,j)==seq_rand_sar(i,j) && seq_evo_sar2(i,j)!=seq_evo_sar(i,j)) || (seq_evo_sar2(i,j)==seq_rand_sar(i,j) && seq_evo_sar1(i,j)!=seq_evo_sar(i,j)))
                        {
                            N_err_All_out_S=N_err_All_out_S+1;
                        }
                        if (seq_evo_sar(i,j)==seq_evo_sar1(i,j))
                        {
                            N_right_SA_out_S=N_right_SA_out_S+1;
                        }
                        else if (seq_evo_sar1(i,j)==seq_rand_sar(i,j))
                        {
                            N_err_SA_out_S=N_err_SA_out_S+1;
                        }
                    }
                }
            }
        }
        SA_out_temp(1,j-1,rand_times)=N_right_SA_out;
        SA_out_temp(2,j-1,rand_times)=N_err_SA_out;
        SA_out_temp(3,j-1,rand_times)=N_right_SA_out_S;
        SA_out_temp(4,j-1,rand_times)=N_err_SA_out_S;
        SA_out_temp(5,j-1,rand_times)=N_right_SA_out-N_right_SA_out_S;
        SA_out_temp(6,j-1,rand_times)=N_err_SA_out-N_err_SA_out_S;
        SA_out_temp(7,j-1,rand_times)=N_uncertainty;
        All_out_temp(1,j-1,rand_times)=N_right_All_out;
        All_out_temp(2,j-1,rand_times)=N_err_All_out;
        All_out_temp(3,j-1,rand_times)=N_right_All_out_S;
        All_out_temp(4,j-1,rand_times)=N_err_All_out_S;
        All_out_temp(5,j-1,rand_times)=N_right_All_out-N_right_All_out_S;
        All_out_temp(6,j-1,rand_times)=N_err_All_out-N_err_All_out_S;
        All_out_temp(7,j-1,rand_times)=N_uncertainty_all;
    }
}
#endif /* confidence_calculation_hpp */
