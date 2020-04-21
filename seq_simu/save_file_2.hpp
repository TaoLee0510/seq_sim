//
//  save_file_2.hpp
//  seq_simu
//
//  Created by Taolee on 4/3/20.
//  Copyright © 2020 Taolee. All rights reserved.
//

#ifndef save_file_2_hpp
#define save_file_2_hpp

#include <stdio.h>
#include <blitz/blitz.h>
#include <blitz/array.h>
#include <cmath>

using namespace blitz;
void save_file_2 (Array<double,3> results0,Array<double,3> results1,Array<double,3> results2,Array<double,3> results3,int D, int j)
{
    char filedir31 [100] = {'\0'};
    sprintf(filedir31, "./results_single_out_%.1d.txt",j);
    FILE * fid31;
    fid31=fopen (filedir31,"w+");
    for(int i=1;i<=D;i++)
    {
        for(int co=1;co<=11;co++)
        {
            if(co<11)
            {
                fprintf(fid31,"%f\t",results0(i,co,j));
            }
            else
            {
                fprintf(fid31,"%f\n",results0(i,co,j));
            }
        }
    }
    fclose(fid31);
    
    
    char filedir5 [100] = {'\0'};
    sprintf(filedir5, "./results_double_out_%.1d.txt",j);
    FILE * fid5;
    fid5=fopen (filedir5,"w+");
    for(int i=1;i<=D;i++)
    {
        for(int co=1;co<=11;co++)
        {
            if(co<11)
            {
                fprintf(fid5,"%f\t",results1(i,co,j));
            }
            else
            {
                fprintf(fid5,"%f\n",results1(i,co,j));
            }
        }
    }
    fclose(fid5);
    
    
    char filedir51 [100] = {'\0'};
    sprintf(filedir51, "./Similarity_%.1d.txt",j);
    FILE * fid51;
    fid51=fopen (filedir51,"w+");
    for(int i=1;i<=D;i++)
    {
        for(int co=1;co<=4;co++)
        {
            if(co<4)
            {
                fprintf(fid51,"%f\t",results2(i,co,j));
            }
            else
            {
                fprintf(fid51,"%f\n",results2(i,co,j));
            }
        }
    }
    fclose(fid51);
    
    char filedir52 [100] = {'\0'};
    sprintf(filedir52, "./Similarity_divergency_%.1d.txt",j);
    FILE * fid52;
    fid52=fopen (filedir52,"w+");
    for(int i=1;i<=D;i++)
    {
        for(int co=1;co<=4;co++)
        {
            if(co<4)
            {
                fprintf(fid52,"%f\t",results3(i,co,j));
            }
            else
            {
                fprintf(fid52,"%f\n",results3(i,co,j));
            }
        }
    }
    fclose(fid52);
}
#endif /* save_file_2_hpp */
