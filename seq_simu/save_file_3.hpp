//
//  save_file_3.hpp
//  seq_simu
//
//  Created by Taolee on 4/13/20.
//  Copyright © 2020 Taolee. All rights reserved.
//

#ifndef save_file_3_hpp
#define save_file_3_hpp

#include <stdio.h>
#include <blitz/blitz.h>
#include <blitz/array.h>
#include <cmath>

using namespace blitz;
void save_file_3 (Array<double,2> results0)
{
        char filedir31 [100] = {'\0'};
        sprintf(filedir31, "./similarity_ref_evo.txt");
        FILE * fid31;
        fid31=fopen (filedir31,"w+");
        for(int co=1;co<=4;co++)
        {
            if(co<4)
            {
                fprintf(fid31,"%f\t",results0(1,co));
            }
            else
            {
                fprintf(fid31,"%f\n",results0(1,co));
            }
        }
        fclose(fid31);
}

#endif /* save_file_3_hpp */
