//
//  accuracy.hpp
//  seq_simu
//
//  Created by Taolee on 4/24/20.
//  Copyright © 2020 Taolee. All rights reserved.
//

#ifndef accuracy_hpp
#define accuracy_hpp

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
void accuracy (Array<int, 2> codon, Array<int, 2> seq_orig_sar, Array<int, 2> seq_evo_sar,Array<int, 2> seq_evo_sar1,Array<int, 2> seq_evo_sar2, Array<int, 2> seq_rand_sar, Array<int, 2> seq_rand_sar1, Array<int, 2> seq_rand_sar0, Array<int, 2> seq_rand_sar01, Array<double, 3> &SA_out_temp, Array<double, 3> &All_out_temp ,Array<double, 3> &Similarity_divergency,Array<double, 3> &Similarity_divergency1, Array<double, 3> &Similarity,Array<double, 3> &Similarity1, int j, int rand_times, int seq_length,int duplic, double NS_site_number, double S_site_number)
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
        int NS_similarity01=0;
        int S_similarity01=0;
        int NS_similarity00=0;
        int S_similarity00=0;
        int NS_similarity1=0;
        int S_similarity1=0;
        int NS_similarity10=0;
        int S_similarity10=0;
        int NS_similarity11=0;
        int S_similarity11=0;
        int NS_similarity3=0;
        int S_similarity3=0;
        int S_uncertainty=0;
        int S_uncertainty_all=0;
        for (int i=1;i<=seq_length;i++)
        {
            int temp=0;
            int temp1=0;
            if (seq_rand_sar(i,j)!=seq_rand_sar1(i,j))
            {
                if (seq_rand_sar0(i,j)==seq_rand_sar(i,j))
                {
                    temp=seq_rand_sar0(i,j);
                }
                if (seq_rand_sar0(i,j)==seq_rand_sar1(i,j))
                {
                    temp=seq_rand_sar0(i,j);
                }
                if (temp==seq_evo_sar(i,j))
                {
                    N_right_SA_out=N_right_SA_out+1;
                }
                if (temp!=seq_evo_sar(i,j))
                {
                    N_err_SA_out=N_err_SA_out+1;
                }
                if (seq_rand_sar0(i,j)!=seq_rand_sar(i,j) && seq_rand_sar0(i,j)!=seq_rand_sar1(i,j))
                {
                    N_uncertainty=N_uncertainty+1;
                }
                if (seq_rand_sar0(i,j)==seq_rand_sar01(i,j))
                {
                    temp1=seq_rand_sar0(i,j);
                    if (temp1==seq_rand_sar(i,j))
                    {
                        temp=temp1;
                    }
                    if (temp1==seq_rand_sar1(i,j))
                    {
                        temp=temp1;
                    }
                    if (temp==seq_evo_sar(i,j))
                    {
                        N_right_All_out=N_right_All_out+1;
                    }
                    else if (temp!=seq_evo_sar(i,j))
                    {
                        N_err_All_out=N_err_All_out+1;
                    }
                    if (temp1!=seq_rand_sar(i,j) && temp1!=seq_rand_sar1(i,j))
                    {
                        N_uncertainty_all=N_uncertainty_all+1;
                    }
                }
                else
                {
                    N_uncertainty_all=N_uncertainty_all+1;
                }
                if (i%3==0)
                {
                    Array<int, 2> seqsar(seq_length,1,FortranArray<2>());
                    seqsar=0;
                    seqsar(all,1)=seq_rand_sar1(all,j);
                    int s=NS_judge(i,seqsar,seq_rand_sar,codon,j);
                    if(s==1)
                    {
                        if (seq_rand_sar0(i,j)==seq_rand_sar(i,j))
                        {
                            temp=seq_rand_sar0(i,j);
                        }
                        if (seq_rand_sar0(i,j)==seq_rand_sar1(i,j))
                        {
                            temp=seq_rand_sar0(i,j);
                        }
                        if (temp==seq_evo_sar(i,j))
                        {
                            N_right_SA_out_S=N_right_SA_out_S+1;
                        }
                        if (temp!=seq_evo_sar(i,j))
                        {
                            N_err_SA_out_S=N_err_SA_out_S+1;
                        }
                        if (seq_rand_sar0(i,j)!=seq_rand_sar(i,j) && seq_rand_sar0(i,j)!=seq_rand_sar1(i,j))
                        {
                            S_uncertainty=S_uncertainty+1;
                        }
                        if (seq_rand_sar0(i,j)==seq_rand_sar01(i,j))
                        {
                            temp1=seq_rand_sar0(i,j);
                            if (temp1==seq_rand_sar(i,j))
                            {
                                temp=temp1;
                            }
                            if (temp1==seq_rand_sar1(i,j))
                            {
                                temp=temp1;
                            }
                            if (temp==seq_evo_sar(i,j))
                            {
                                N_right_All_out_S=N_right_All_out_S+1;
                            }
                            else if (temp!=seq_evo_sar(i,j))
                            {
                                N_err_All_out_S=N_err_All_out_S+1;
                            }
                            if (temp1!=seq_rand_sar(i,j) && temp1!=seq_rand_sar1(i,j))
                            {
                                S_uncertainty_all=S_uncertainty_all+1;
                            }
                        }
                    }
                }
            }
            if (seq_rand_sar0(i,j)!=seq_rand_sar(i,j))
            {
                if (i%3==0)
                {
                    Array<int, 2> seqsar(seq_length,1,FortranArray<2>());
                    seqsar=0;
                    seqsar(all,1)=seq_rand_sar0(all,j);
                    int s=NS_judge(i,seqsar,seq_rand_sar,codon,j);
                    if (s==1)
                    {
                        S_similarity00=S_similarity00+1;
                    }
                    else
                    {
                        NS_similarity00=NS_similarity00+1;
                    }
                }
                else
                {
                    NS_similarity00=NS_similarity00+1;
                }
            }
            if (seq_rand_sar0(i,j)!=seq_rand_sar1(i,j))
            {
                if (i%3==0)
                {
                    Array<int, 2> seqsar(seq_length,1,FortranArray<2>());
                    seqsar=0;
                    seqsar(all,1)=seq_rand_sar0(all,j);
                    int s=NS_judge(i,seqsar,seq_rand_sar1,codon,j);
                    if (s==1)
                    {
                        S_similarity01=S_similarity01+1;
                    }
                    else
                    {
                        NS_similarity01=NS_similarity01+1;
                    }
                }
                else
                {
                    NS_similarity01=NS_similarity01+1;
                }
            }
            if (seq_rand_sar01(i,j)!=seq_rand_sar(i,j))
            {
                if (i%3==0)
                {
                    Array<int, 2> seqsar(seq_length,1,FortranArray<2>());
                    seqsar=0;
                    seqsar(all,1)=seq_rand_sar01(all,j);
                    int s=NS_judge(i,seqsar,seq_rand_sar,codon,j);
                    if (s==1)
                    {
                        S_similarity10=S_similarity10+1;
                    }
                    else
                    {
                        NS_similarity10=NS_similarity10+1;
                    }
                }
                else
                {
                    NS_similarity10=NS_similarity10+1;
                }
            }
            if (seq_rand_sar01(i,j)!=seq_rand_sar1(i,j))
            {
                if (i%3==0)
                {
                    Array<int, 2> seqsar(seq_length,1,FortranArray<2>());
                    seqsar=0;
                    seqsar(all,1)=seq_rand_sar01(all,j);
                    int s=NS_judge(i,seqsar,seq_rand_sar1,codon,j);
                    if (s==1)
                    {
                        S_similarity11=S_similarity11+1;
                    }
                    else
                    {
                        NS_similarity11=NS_similarity11+1;
                    }
                }
                else
                {
                    NS_similarity11=NS_similarity11+1;
                }
            }
            if (seq_rand_sar(i,j)!=seq_rand_sar1(i,j))
            {
                if (i%3==0)
                {
                    Array<int, 2> seqsar(seq_length,1,FortranArray<2>());
                    seqsar=0;
                    seqsar(all,1)=seq_rand_sar(all,j);
                    int s=NS_judge(i,seqsar,seq_rand_sar1,codon,j);
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
            if (seq_orig_sar(i,1)!=seq_evo_sar(i,j))
            {
                if (i%3==0)
                {
                    int s=NS_judge(i,seq_orig_sar,seq_evo_sar,codon,j);
                    if (s==1)
                    {
                        S_similarity3=S_similarity3+1;
                    }
                    else
                    {
                        NS_similarity3=NS_similarity3+1;
                    }
                }
                else
                {
                    NS_similarity3=NS_similarity3+1;
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
        SA_out_temp(8,j-1,rand_times)=S_uncertainty;
        SA_out_temp(9,j-1,rand_times)=N_uncertainty-S_uncertainty;
        All_out_temp(1,j-1,rand_times)=N_right_All_out;
        All_out_temp(2,j-1,rand_times)=N_err_All_out;
        All_out_temp(3,j-1,rand_times)=N_right_All_out_S;
        All_out_temp(4,j-1,rand_times)=N_err_All_out_S;
        All_out_temp(5,j-1,rand_times)=N_right_All_out-N_right_All_out_S;
        All_out_temp(6,j-1,rand_times)=N_err_All_out-N_err_All_out_S;
        All_out_temp(7,j-1,rand_times)=N_uncertainty_all;
        All_out_temp(8,j-1,rand_times)=S_uncertainty_all;
        All_out_temp(9,j-1,rand_times)=N_uncertainty_all-S_uncertainty_all;
        Similarity(1,j-1,duplic)=j;
        double sum5=S_similarity00+NS_similarity00;
        double sum51=S_similarity01+NS_similarity01;
        Similarity(2,j-1,duplic)=(sum5+sum51)/(2*seq_length);
        Similarity(3,j-1,duplic)=(S_similarity00+S_similarity01)/(2*S_site_number);
        Similarity(4,j-1,duplic)=(NS_similarity00+NS_similarity01)/(2*NS_site_number);
        Similarity_divergency(1,j-1,duplic)=j;
        double sum6=S_similarity1+NS_similarity1;
        Similarity_divergency(2,j-1,duplic)=sum6/seq_length;
        Similarity_divergency(3,j-1,duplic)=S_similarity1/S_site_number;
        Similarity_divergency(4,j-1,duplic)=NS_similarity1/NS_site_number;
        Similarity_divergency1(1,j-1,duplic)=j;
        double sum7=S_similarity10+NS_similarity10;
        double sum71=S_similarity11+NS_similarity11;
        Similarity_divergency1(2,j-1,duplic)=(sum7+sum71)/(2*seq_length);
        Similarity_divergency1(3,j-1,duplic)=(S_similarity10+S_similarity11)/(2*S_site_number);
        Similarity_divergency1(4,j-1,duplic)=(NS_similarity10+NS_similarity11)/(2*NS_site_number);
        Similarity1(1,j-1,duplic)=j;
        double sum8=S_similarity3+NS_similarity3;
        Similarity1(2,j-1,duplic)=sum8/seq_length;
        Similarity1(3,j-1,duplic)=S_similarity3/S_site_number;
        Similarity1(4,j-1,duplic)=NS_similarity3/NS_site_number;
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
        int S_uncertainty=0;
        int S_uncertainty_all=0;
        for (int i=1;i<=seq_length;i++)
        {
           int temp=0;
            int temp1=0;
            if (seq_rand_sar(i,j)!=seq_rand_sar1(i,j))
            {
                if (seq_rand_sar0(i,j)==seq_rand_sar(i,j))
                {
                    temp=seq_rand_sar0(i,j);
                }
                if (seq_rand_sar0(i,j)==seq_rand_sar1(i,j))
                {
                    temp=seq_rand_sar0(i,j);
                }
                if (temp==seq_evo_sar(i,j))
                {
                    N_right_SA_out=N_right_SA_out+1;
                }
                if (temp!=seq_evo_sar(i,j))
                {
                    N_err_SA_out=N_err_SA_out+1;
                }
                if (seq_rand_sar0(i,j)!=seq_rand_sar(i,j) && seq_rand_sar0(i,j)!=seq_rand_sar1(i,j))
                {
                    N_uncertainty=N_uncertainty+1;
                }
                if (seq_rand_sar0(i,j)==seq_rand_sar01(i,j))
                {
                    temp1=seq_rand_sar0(i,j);
                    if (temp1==seq_rand_sar(i,j))
                    {
                        temp=temp1;
                    }
                    if (temp1==seq_rand_sar1(i,j))
                    {
                        temp=temp1;
                    }
                    if (temp==seq_evo_sar(i,j))
                    {
                        N_right_All_out=N_right_All_out+1;
                    }
                    else if (temp!=seq_evo_sar(i,j))
                    {
                        N_err_All_out=N_err_All_out+1;
                    }
                    if (temp1!=seq_rand_sar(i,j) && temp1!=seq_rand_sar1(i,j))
                    {
                        N_uncertainty_all=N_uncertainty_all+1;
                    }
                }
                else
                {
                    N_uncertainty_all=N_uncertainty_all+1;
                }
                if (i%3==0)
                {
                    Array<int, 2> seqsar(seq_length,1,FortranArray<2>());
                    seqsar=0;
                    seqsar(all,1)=seq_rand_sar1(all,j);
                    int s=NS_judge(i,seqsar,seq_rand_sar,codon,j);
                    if(s==1)
                    {
                        if (seq_rand_sar0(i,j)==seq_rand_sar(i,j))
                        {
                            temp=seq_rand_sar0(i,j);
                        }
                        if (seq_rand_sar0(i,j)==seq_rand_sar1(i,j))
                        {
                            temp=seq_rand_sar0(i,j);
                        }
                        if (temp==seq_evo_sar(i,j))
                        {
                            N_right_SA_out_S=N_right_SA_out_S+1;
                        }
                        else if (temp!=seq_evo_sar(i,j))
                        {
                            N_err_SA_out_S=N_err_SA_out_S+1;
                        }
                        if (seq_rand_sar0(i,j)!=seq_rand_sar(i,j) && seq_rand_sar0(i,j)!=seq_rand_sar1(i,j))
                        {
                            S_uncertainty=S_uncertainty+1;
                        }
                        if (seq_rand_sar0(i,j)==seq_rand_sar01(i,j))
                        {
                            temp1=seq_rand_sar0(i,j);
                            if (temp1==seq_rand_sar(i,j))
                            {
                                temp=temp1;
                            }
                            if (temp1==seq_rand_sar1(i,j))
                            {
                                temp=temp1;
                            }
                            if (temp==seq_evo_sar(i,j))
                            {
                                N_right_All_out_S=N_right_All_out_S+1;
                            }
                            else if (temp!=seq_evo_sar(i,j))
                            {
                                N_err_All_out_S=N_err_All_out_S+1;
                            }
                            if (temp1!=seq_rand_sar(i,j) && temp1!=seq_rand_sar1(i,j))
                            {
                                S_uncertainty_all=S_uncertainty_all+1;
                            }
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
        SA_out_temp(8,j-1,rand_times)=S_uncertainty;
        SA_out_temp(9,j-1,rand_times)=N_uncertainty-S_uncertainty;
        All_out_temp(1,j-1,rand_times)=N_right_All_out;
        All_out_temp(2,j-1,rand_times)=N_err_All_out;
        All_out_temp(3,j-1,rand_times)=N_right_All_out_S;
        All_out_temp(4,j-1,rand_times)=N_err_All_out_S;
        All_out_temp(5,j-1,rand_times)=N_right_All_out-N_right_All_out_S;
        All_out_temp(6,j-1,rand_times)=N_err_All_out-N_err_All_out_S;
        All_out_temp(7,j-1,rand_times)=N_uncertainty_all;
        All_out_temp(8,j-1,rand_times)=S_uncertainty_all;
        All_out_temp(9,j-1,rand_times)=N_uncertainty_all-S_uncertainty_all;
    }
}

#endif /* accuracy_hpp */

