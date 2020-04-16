//
//  extra_results.hpp
//  seq_simu
//
//  Created by Taolee on 3/29/20.
//  Copyright © 2020 Taolee. All rights reserved.
//

#ifndef extra_results_hpp
#define extra_results_hpp

#include <stdio.h>
#include <blitz/blitz.h>
#include <blitz/array.h>
#include <math.h>

using namespace std;
using namespace blitz;
void extral_results(Array<double, 3> SA_out_temp, Array<double, 3> All_out_temp,Array<double, 3> &SA_out,Array<double, 3> &All_out,int D, int duplic)
{
    for (int j=1;j<=D;j++)
    {
        Range all = Range::all();
        SA_out(1,j,duplic)=mean(SA_out_temp(1,j,all));
        SA_out(2,j,duplic)=mean(SA_out_temp(2,j,all));
        double sum1=SA_out(1,j,duplic)+SA_out(2,j,duplic);
        if (sum1==0)
        {
            SA_out(3,j,duplic)=0;
        }
        else
        {
            SA_out(3,j,duplic)=SA_out(2,j,duplic)/sum1;
        }
        SA_out(4,j,duplic)=mean(SA_out_temp(3,j,all));
        SA_out(5,j,duplic)=mean(SA_out_temp(4,j,all));
        double sum2=SA_out(4,j,duplic)+SA_out(5,j,duplic);
        if (sum2==0)
        {
            SA_out(6,j,duplic)=0;
        }
        else
        {
            SA_out(6,j,duplic)=SA_out(5,j,duplic)/sum2;
        }
        SA_out(7,j,duplic)=mean(SA_out_temp(5,j,all));
        SA_out(8,j,duplic)=mean(SA_out_temp(6,j,all));
        double sum3=SA_out(7,j,duplic)+SA_out(8,j,duplic);
        if (sum3==0)
        {
            SA_out(9,j,duplic)=0;
        }
        else
        {
            SA_out(9,j,duplic)=SA_out(8,j,duplic)/sum3;
        }
        SA_out(10,j,duplic)=mean(SA_out_temp(7,j,all));
        double sum4=SA_out(10,j,duplic)+SA_out(1,j,duplic)+SA_out(2,j,duplic);
        if (sum4==0)
        {
            SA_out(11,j,duplic)=0;
        }
        else
        {
            SA_out(11,j,duplic)=SA_out(10,j,duplic)/sum4;
        }
        All_out(1,j,duplic)=mean(All_out_temp(1,j,all));
        All_out(2,j,duplic)=mean(All_out_temp(2,j,all));
        sum1=All_out(1,j,duplic)+All_out(2,j,duplic);
        if (sum1==0)
        {
            All_out(3,j,duplic)=0;
        }
        else
        {
            All_out(3,j,duplic)=All_out(2,j,duplic)/sum1;
        }
        All_out(4,j,duplic)=mean(All_out_temp(3,j,all));
        All_out(5,j,duplic)=mean(All_out_temp(4,j,all));
        sum2=All_out(4,j,duplic)+All_out(5,j,duplic);
        if (sum2==0)
        {
            All_out(6,j,duplic)=0;
        }
        else
        {
            All_out(6,j,duplic)=All_out(5,j,duplic)/sum2;
        }
        All_out(7,j,duplic)=mean(All_out_temp(5,j,all));
        All_out(8,j,duplic)=mean(All_out_temp(6,j,all));
        sum3=All_out(7,j,duplic)+All_out(8,j,duplic);
        if (sum3==0)
        {
            All_out(9,j,duplic)=0;
        }
        else
        {
            All_out(9,j,duplic)=All_out(8,j,duplic)/sum3;
        }
        All_out(10,j,duplic)=mean(All_out_temp(7,j,all));
        sum4=All_out(10,j,duplic)+All_out(j,1,duplic)+All_out(2,j,duplic);
        if (sum4==0)
        {
            All_out(11,j,duplic)=0;
        }
        else
        {
            All_out(11,j,duplic)=All_out(10,j,duplic)/sum4;
        }
    }
}
#endif /* extra_results_hpp */
