//
//  confidence_calculation.hpp
//  seq_simu
//
//  Created by Taolee on 4/16/20.
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

void confidence_calculation (Array<int, 2> codon, Array<int, 2> seq_orig_sar, Array<int, 2> seq_evo_sar,Array<int, 2> seq_evo_sar1,Array<int, 2> seq_evo_sar2, Array<int, 2> seq_rand_sar,Array<double, 3> &SA_out_temp, Array<double, 3> &All_out_temp ,Array<double, 3> &Similarity_divergency,Array<double, 3> &Similarity,int j, int rand_times, int seq_length,int duplic)
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
            if (seq_rand_sar(j,i)!=seq_evo_sar(j,i))
            {
                if ((seq_evo_sar1(j,i)==seq_evo_sar2(j,i) && seq_evo_sar1(j,i)==seq_evo_sar(j,i)) || (seq_evo_sar1(j,i)==seq_evo_sar(j,i) && seq_evo_sar2(j,i)!=seq_rand_sar(j,i))|| (seq_evo_sar2(j,i)==seq_evo_sar(j,i) && seq_evo_sar1(j,i)!=seq_rand_sar(j,i)))
                {
                    N_right_All_out=N_right_All_out+1;
                }
                else if ((seq_evo_sar1(j,i)==seq_evo_sar2(j,i) && seq_evo_sar1(j,i)==seq_rand_sar(j,i)) || (seq_evo_sar1(j,i)==seq_rand_sar(j,i) && seq_evo_sar2(j,i)!=seq_evo_sar(j,i)) || (seq_evo_sar2(j,i)==seq_rand_sar(j,i) && seq_evo_sar1(j,i)!=seq_evo_sar(j,i)))
                {
                    N_err_All_out=N_err_All_out+1;
                }
                else if ((seq_evo_sar2(j,i)!=seq_rand_sar(j,i) && seq_evo_sar2(j,i)!=seq_evo_sar(j,i)) && (seq_evo_sar1(j,i)!=seq_rand_sar(j,i) && seq_evo_sar1(j,i)!=seq_evo_sar(j,i)))
                {
                    N_uncertainty_all=N_uncertainty_all+1;
                }
                if (seq_evo_sar(j,i)==seq_evo_sar1(j,i))
                {
                    N_right_SA_out=N_right_SA_out+1;
                }
                else if (seq_evo_sar1(j,i)==seq_rand_sar(j,i))
                {
                    N_err_SA_out=N_err_SA_out+1;
                }
                else
                {
                    N_uncertainty=N_uncertainty+1;
                }
                total_mutations=total_mutations+1;
            }
            if (i%3==0)
            {
                if (seq_rand_sar(j,i)!=seq_evo_sar(j,i))
                {
                    Array<int, 2> seqsar(1,seq_length,FortranArray<2>());//seq_rand
                    seqsar=0;
                    seqsar(1,all)=seq_evo_sar(j,all);
                    int s=NS_judge(i,seqsar,seq_rand_sar,codon,j);
                    if(s==1)
                    {
                        if ((seq_evo_sar1(j,i)==seq_evo_sar2(j,i) && seq_evo_sar1(j,i)==seq_evo_sar(j,i)) || (seq_evo_sar1(j,i)==seq_evo_sar(j,i) && seq_evo_sar2(j,i)!=seq_rand_sar(j,i))|| (seq_evo_sar2(j,i)==seq_evo_sar(j,i) && seq_evo_sar1(j,i)!=seq_rand_sar(j,i)))
                        {
                            N_right_All_out_S=N_right_All_out_S+1;
                        }
                        else if ((seq_evo_sar1(j,i)==seq_evo_sar2(j,i) && seq_evo_sar1(j,i)==seq_rand_sar(j,i)) || (seq_evo_sar1(j,i)==seq_rand_sar(j,i) && seq_evo_sar2(j,i)!=seq_evo_sar(j,i)) || (seq_evo_sar2(j,i)==seq_rand_sar(j,i) && seq_evo_sar1(j,i)!=seq_evo_sar(j,i)))
                        {
                            N_err_All_out_S=N_err_All_out_S+1;
                        }
                        if (seq_evo_sar(j,i)==seq_evo_sar1(j,i))
                        {
                            N_right_SA_out_S=N_right_SA_out_S+1;
                        }
                        else if (seq_evo_sar1(j,i)==seq_rand_sar(j,i))
                        {
                            N_err_SA_out_S=N_err_SA_out_S+1;
                        }
                    }
                }
            }
            if (seq_evo_sar1(j,i)!=seq_evo_sar(j,i))
            {
                if (i%3==0)
                {
                    Array<int, 2> seqsar(1,seq_length,FortranArray<2>());//seq_rand
                    seqsar=0;
                    seqsar(1,all)=seq_evo_sar1(j,all);
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
            if (seq_rand_sar(j,i)!=seq_evo_sar(j,i))
            {
                if (i%3==0)
                {
                    Array<int, 2> seqsar(1,seq_length,FortranArray<2>());//seq_rand
                    seqsar=0;
                    seqsar(1,all)=seq_rand_sar(j,all);
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
        
        SA_out_temp(j-1,1,rand_times)=N_right_SA_out;
        SA_out_temp(j-1,2,rand_times)=N_err_SA_out;
        SA_out_temp(j-1,3,rand_times)=N_right_SA_out_S;
        SA_out_temp(j-1,4,rand_times)=N_err_SA_out_S;
        SA_out_temp(j-1,5,rand_times)=N_right_SA_out-N_right_SA_out_S;
        SA_out_temp(j-1,6,rand_times)=N_err_SA_out-N_err_SA_out_S;
        SA_out_temp(j-1,7,rand_times)=N_uncertainty;
        
        All_out_temp(j-1,1,rand_times)=N_right_All_out;
        All_out_temp(j-1,2,rand_times)=N_err_All_out;
        All_out_temp(j-1,3,rand_times)=N_right_All_out_S;
        All_out_temp(j-1,4,rand_times)=N_err_All_out_S;
        All_out_temp(j-1,5,rand_times)=N_right_All_out-N_right_All_out_S;
        All_out_temp(j-1,6,rand_times)=N_err_All_out-N_err_All_out_S;
        All_out_temp(j-1,7,rand_times)=N_uncertainty_all;
        
        cout << "Runing times: "<<duplic<< "     Confidence run times:    "<< rand_times <<"     Confidence run: "<< j <<endl;
        
        Similarity(j-1,1,duplic)=j;
        double NS_site_number=seq_length*3.5/4.5;
        double S_site_number=seq_length/4.5;
        double sum5=S_similarity+NS_similarity;
        Similarity(j-1,2,duplic)=sum5/seq_length;
        Similarity(j-1,3,duplic)=S_similarity/S_site_number;
        Similarity(j-1,4,duplic)=NS_similarity/NS_site_number;
        Similarity_divergency(j-1,1,duplic)=j;
        double NS_site_number1=seq_length*3.5/4.5;
        double S_site_number1=seq_length/4.5;
        double sum6=S_similarity1+NS_similarity1;
        Similarity_divergency(j-1,2,duplic)=sum6/seq_length;
        Similarity_divergency(j-1,3,duplic)=S_similarity1/S_site_number1;
        Similarity_divergency(j-1,4,duplic)=NS_similarity1/NS_site_number1;
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
            if (seq_rand_sar(j,i)!=seq_evo_sar(j,i))
            {
                if ((seq_evo_sar1(j,i)==seq_evo_sar2(j,i) && seq_evo_sar1(j,i)==seq_evo_sar(j,i)) || (seq_evo_sar1(j,i)==seq_evo_sar(j,i) && seq_evo_sar2(j,i)!=seq_rand_sar(j,i))|| (seq_evo_sar2(j,i)==seq_evo_sar(j,i) && seq_evo_sar1(j,i)!=seq_rand_sar(j,i)))
                {
                    N_right_All_out=N_right_All_out+1;
                }
                else if ((seq_evo_sar1(j,i)==seq_evo_sar2(j,i) && seq_evo_sar1(j,i)==seq_rand_sar(j,i)) || (seq_evo_sar1(j,i)==seq_rand_sar(j,i) && seq_evo_sar2(j,i)!=seq_evo_sar(j,i)) || (seq_evo_sar2(j,i)==seq_rand_sar(j,i) && seq_evo_sar1(j,i)!=seq_evo_sar(j,i)))
                {
                    N_err_All_out=N_err_All_out+1;
                }
                else if ((seq_evo_sar2(j,i)!=seq_rand_sar(j,i) && seq_evo_sar2(j,i)!=seq_evo_sar(j,i)) && (seq_evo_sar1(j,i)!=seq_rand_sar(j,i) && seq_evo_sar1(j,i)!=seq_evo_sar(j,i)))
                {
                    N_uncertainty_all=N_uncertainty_all+1;
                }
                if (seq_evo_sar(j,i)==seq_evo_sar1(j,i))
                {
                    N_right_SA_out=N_right_SA_out+1;
                }
                else if (seq_evo_sar1(j,i)==seq_rand_sar(j,i))
                {
                    N_err_SA_out=N_err_SA_out+1;
                }
                else
                {
                    N_uncertainty=N_uncertainty+1;
                }
                total_mutations=total_mutations+1;
            }
            if (i%3==0)
            {
                if (seq_rand_sar(j,i)!=seq_evo_sar(j,i))
                {
                    Array<int, 2> seqsar(1,seq_length,FortranArray<2>());//seq_rand
                    seqsar=0;
                    seqsar(1,all)=seq_evo_sar(j,all);
                    int s=NS_judge(i,seqsar,seq_rand_sar,codon,j);
                    if(s==1)
                    {
                        if ((seq_evo_sar1(j,i)==seq_evo_sar2(j,i) && seq_evo_sar1(j,i)==seq_evo_sar(j,i)) || (seq_evo_sar1(j,i)==seq_evo_sar(j,i) && seq_evo_sar2(j,i)!=seq_rand_sar(j,i))|| (seq_evo_sar2(j,i)==seq_evo_sar(j,i) && seq_evo_sar1(j,i)!=seq_rand_sar(j,i)))
                        {
                            N_right_All_out_S=N_right_All_out_S+1;
                        }
                        else if ((seq_evo_sar1(j,i)==seq_evo_sar2(j,i) && seq_evo_sar1(j,i)==seq_rand_sar(j,i)) || (seq_evo_sar1(j,i)==seq_rand_sar(j,i) && seq_evo_sar2(j,i)!=seq_evo_sar(j,i)) || (seq_evo_sar2(j,i)==seq_rand_sar(j,i) && seq_evo_sar1(j,i)!=seq_evo_sar(j,i)))
                        {
                            N_err_All_out_S=N_err_All_out_S+1;
                        }
                        if (seq_evo_sar(j,i)==seq_evo_sar1(j,i))
                        {
                            N_right_SA_out_S=N_right_SA_out_S+1;
                        }
                        else if (seq_evo_sar1(j,i)==seq_rand_sar(j,i))
                        {
                            N_err_SA_out_S=N_err_SA_out_S+1;
                        }
                    }
                }
            }
        }
        SA_out_temp(j-1,1,rand_times)=N_right_SA_out;
        SA_out_temp(j-1,2,rand_times)=N_err_SA_out;
        SA_out_temp(j-1,3,rand_times)=N_right_SA_out_S;
        SA_out_temp(j-1,4,rand_times)=N_err_SA_out_S;
        SA_out_temp(j-1,5,rand_times)=N_right_SA_out-N_right_SA_out_S;
        SA_out_temp(j-1,6,rand_times)=N_err_SA_out-N_err_SA_out_S;
        SA_out_temp(j-1,7,rand_times)=N_uncertainty;
        
        All_out_temp(j-1,1,rand_times)=N_right_All_out;
        All_out_temp(j-1,2,rand_times)=N_err_All_out;
        All_out_temp(j-1,3,rand_times)=N_right_All_out_S;
        All_out_temp(j-1,4,rand_times)=N_err_All_out_S;
        All_out_temp(j-1,5,rand_times)=N_right_All_out-N_right_All_out_S;
        All_out_temp(j-1,6,rand_times)=N_err_All_out-N_err_All_out_S;
        All_out_temp(j-1,7,rand_times)=N_uncertainty_all;
        
        cout << "Runing times: "<<duplic<< "     Confidence run times:    "<< rand_times <<"     Confidence run: "<< j <<endl;
    }
}

#endif /* confidence_calculation_hpp */
