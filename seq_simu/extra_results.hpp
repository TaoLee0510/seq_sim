//
//  extra_results.hpp
//  seq_simu
//
//  Created by Taolee on 4/15/20.
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
        SA_out(j,1,duplic)=mean(SA_out_temp(j,1,all));
        SA_out(j,2,duplic)=mean(SA_out_temp(j,2,all));
        double sum1=SA_out(j,1,duplic)+SA_out(j,2,duplic);
        if (sum1==0)
        {
            SA_out(j,3,duplic)=0;
        }
        else
        {
            SA_out(j,3,duplic)=SA_out(j,2,duplic)/sum1;
        }
        SA_out(j,4,duplic)=mean(SA_out_temp(j,3,all));
        SA_out(j,5,duplic)=mean(SA_out_temp(j,4,all));
        double sum2=SA_out(j,4,duplic)+SA_out(j,5,duplic);
        if (sum2==0)
        {
            SA_out(j,6,duplic)=0;
        }
        else
        {
            SA_out(j,6,duplic)=SA_out(j,5,duplic)/sum2;
        }
        SA_out(j,7,duplic)=mean(SA_out_temp(j,5,all));
        SA_out(j,8,duplic)=mean(SA_out_temp(j,6,all));
        double sum3=SA_out(j,7,duplic)+SA_out(j,8,duplic);
        if (sum3==0)
        {
            SA_out(j,9,duplic)=0;
        }
        else
        {
            SA_out(j,9,duplic)=SA_out(j,8,duplic)/sum3;
        }
        SA_out(j,10,duplic)=mean(SA_out_temp(j,7,all));
        double sum4=SA_out(j,10,duplic)+SA_out(j,1,duplic)+SA_out(j,2,duplic);
        if (sum4==0)
        {
            SA_out(j,11,duplic)=0;
        }
        else
        {
            SA_out(j,11,duplic)=SA_out(j,10,duplic)/sum4;
        }
        
        
        All_out(j,1,duplic)=mean(All_out_temp(j,1,all));
        All_out(j,2,duplic)=mean(All_out_temp(j,2,all));
        sum1=All_out(j,1,duplic)+All_out(j,2,duplic);
        if (sum1==0)
        {
            All_out(j,3,duplic)=0;
        }
        else
        {
            All_out(j,3,duplic)=All_out(j,2,duplic)/sum1;
        }
        All_out(j,4,duplic)=mean(All_out_temp(j,3,all));
        All_out(j,5,duplic)=mean(All_out_temp(j,4,all));
        sum2=All_out(j,4,duplic)+All_out(j,5,duplic);
        if (sum2==0)
        {
            All_out(j,6,duplic)=0;
        }
        else
        {
            All_out(j,6,duplic)=All_out(j,5,duplic)/sum2;
        }
        All_out(j,7,duplic)=mean(All_out_temp(j,5,all));
        All_out(j,8,duplic)=mean(All_out_temp(j,6,all));
        sum3=All_out(j,7,duplic)+All_out(j,8,duplic);
        if (sum3==0)
        {
            All_out(j,9,duplic)=0;
        }
        else
        {
            All_out(j,9,duplic)=All_out(j,8,duplic)/sum3;
        }
        All_out(j,10,duplic)=mean(All_out_temp(j,7,all));
        sum4=All_out(j,10,duplic)+All_out(j,1,duplic)+All_out(j,2,duplic);
        if (sum4==0)
        {
            All_out(j,11,duplic)=0;
        }
        else
        {
            All_out(j,11,duplic)=All_out(j,10,duplic)/sum4;
        }
    }
}


#endif /* extra_results_hpp */
