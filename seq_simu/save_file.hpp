//
//  save_file.hpp
//  seq_simu
//
//  Created by Taolee on 3/30/20.
//  Copyright © 2020 Taolee. All rights reserved.
//

#ifndef save_file_hpp
#define save_file_hpp
#include <stdio.h>
#include <blitz/blitz.h>
#include <blitz/array.h>
#include <cmath>

using namespace blitz;
void save_file (Array<int,2> seq_1,Array<int,2> seq_11,Array<int,2> seq_12,int duplicates)
{
    char filedir2 [100] = {'\0'};
    sprintf(filedir2, "./Seq_evo_%.1d.txt",duplicates);
    FILE * fid2;
    fid2=fopen (filedir2,"w+");
    int C1 = seq_1.cols();
    int C2 = seq_1.rows();
    for(int cr=1;cr<=C1;cr++)
    {
        for (int cl=1;cl<=C2;cl++)
        {
            
            if(cl<C2)
            {
                if (seq_1(cl,cr)==1)
                {
                    fprintf(fid2,"%s","A");
                }
                else if (seq_1(cl,cr)==2)
                {
                    fprintf(fid2,"%s","T");
                }
                else if (seq_1(cl,cr)==3)
                {
                    fprintf(fid2,"%s","G");
                }
                else if (seq_1(cl,cr)==4)
                {
                    fprintf(fid2,"%s","C");
                }
                else if (seq_1(cl,cr)==0)
                {
                    fprintf(fid2,"%s","-");
                }
            }
            else
            {
                if (seq_1(cl,cr)==1)
                {
                    fprintf(fid2,"%s\n","A");
                }
                else if (seq_1(cl,cr)==2)
                {
                    fprintf(fid2,"%s\n","T");
                }
                else if (seq_1(cl,cr)==3)
                {
                    fprintf(fid2,"%s\n","G");
                }
                else if (seq_1(cl,cr)==4)
                {
                    fprintf(fid2,"%s\n","C");
                }
                else if (seq_1(cl,cr)==0)
                {
                    fprintf(fid2,"%s\n","-");
                }
            }
        }
    }
    fclose(fid2);
    char filedir21 [100] = {'\0'};
    sprintf(filedir21, "./Seq_evo1_%.1d.txt",duplicates);
    FILE * fid21;
    fid21=fopen (filedir21,"w+");
    C1 = seq_11.cols();
    C2 = seq_11.rows();
    for(int cr=1;cr<=C1;cr++)
    {
        for (int cl=1;cl<=C2;cl++)
        {
            
            if(cl<C2)
            {
                if (seq_11(cl,cr)==1)
                {
                    fprintf(fid21,"%s","A");
                }
                else if (seq_11(cl,cr)==2)
                {
                    fprintf(fid21,"%s","T");
                }
                else if (seq_11(cl,cr)==3)
                {
                    fprintf(fid21,"%s","G");
                }
                else if (seq_11(cl,cr)==4)
                {
                    fprintf(fid21,"%s","C");
                }
                else if (seq_11(cl,cr)==0)
                {
                    fprintf(fid21,"%s","-");
                }
            }
            else
            {
                if (seq_11(cl,cr)==1)
                {
                    fprintf(fid21,"%s\n","A");
                }
                else if (seq_11(cl,cr)==2)
                {
                    fprintf(fid21,"%s\n","T");
                }
                else if (seq_11(cl,cr)==3)
                {
                    fprintf(fid21,"%s\n","G");
                }
                else if (seq_11(cl,cr)==4)
                {
                    fprintf(fid21,"%s\n","C");
                }
                else if (seq_11(cl,cr)==0)
                {
                    fprintf(fid21,"%s\n","-");
                }
            }
        }
    }
    fclose(fid21);
    char filedir22 [100] = {'\0'};
    sprintf(filedir22, "./Seq_evo2_%.1d.txt",duplicates);
    FILE * fid22;
    fid22=fopen (filedir22,"w+");
    C1 = seq_12.cols();
    C2 = seq_12.rows();
    for(int cr=1;cr<=C1;cr++)
    {
        for (int cl=1;cl<=C2;cl++)
        {
            
            if(cl<C2)
            {
                if (seq_12(cl,cr)==1)
                {
                    fprintf(fid22,"%s","A");
                }
                else if (seq_12(cl,cr)==2)
                {
                    fprintf(fid22,"%s","T");
                }
                else if (seq_12(cl,cr)==3)
                {
                    fprintf(fid22,"%s","G");
                }
                else if (seq_12(cl,cr)==4)
                {
                    fprintf(fid22,"%s","C");
                }
                else if (seq_12(cl,cr)==0)
                {
                    fprintf(fid22,"%s","-");
                }
            }
            else
            {
                if (seq_12(cl,cr)==1)
                {
                    fprintf(fid22,"%s\n","A");
                }
                else if (seq_12(cl,cr)==2)
                {
                    fprintf(fid22,"%s\n","T");
                }
                else if (seq_12(cl,cr)==3)
                {
                    fprintf(fid22,"%s\n","G");
                }
                else if (seq_12(cl,cr)==4)
                {
                    fprintf(fid22,"%s\n","C");
                }
                else if (seq_12(cl,cr)==0)
                {
                    fprintf(fid22,"%s\n","-");
                }
            }
        }
    }
    fclose(fid22);
}
#endif /* save_file_hpp */
