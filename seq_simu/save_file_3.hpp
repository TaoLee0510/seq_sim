//
//  save_file_3.hpp
//  seq_simu
//
//  Created by Taolee on 3/29/20.
//  Copyright © 2020 Taolee. All rights reserved.
//

#ifndef save_file_3_hpp
#define save_file_3_hpp

#include <stdio.h>
#include <blitz/blitz.h>
#include <blitz/array.h>
#include <cmath>
using namespace blitz;
void save_file_3 (Array<double,2> results0, Array<double,2> results1,int times)
{
    char filedir31 [100] = {'\0'};
    sprintf(filedir31, "./similarity_ref_evo.txt");
    FILE * fid31;
    fid31=fopen (filedir31,"w+");
    int ii=times+1;
    for (int i=1;i<=ii;i++)
    {
        for(int co=1;co<=4;co++)
        {
            if(co<4)
            {
                fprintf(fid31,"%f\t",results0(i,co));
            }
            else
            {
                fprintf(fid31,"%f\n",results0(i,co));
            }
        }
    }
    fclose(fid31);
    
    char filedir3 [100] = {'\0'};
    sprintf(filedir3, "./mutation.txt");
    FILE * fid3;
    fid3=fopen (filedir3,"w+");
    for (int ro=1;ro<=times;ro++)
    {
        for(int co=1;co<=10;co++)
        {
            if(co<10)
            {
                fprintf(fid3,"%f\t",results1(ro,co));
            }
            else
            {
                fprintf(fid3,"%f\n",results1(ro,co));
            }
        }
    }
    fclose(fid3);
}
#endif /* save_file_3_hpp */
