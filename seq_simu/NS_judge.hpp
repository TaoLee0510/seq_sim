//
//  NS_judge.hpp
//  seq_simu
//
//  Created by Taolee on 4/1/20.
//  Copyright © 2020 Taolee. All rights reserved.
//

#ifndef NS_judge_hpp
#define NS_judge_hpp

#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <ctime>
#include <blitz/blitz.h>
#include <blitz/array.h>
using namespace std;
using namespace blitz;

int NS_judge(int i,Array<int, 2> seq1,Array<int, 2> seq2,Array<int, 2> codon,int DAY)
{
    Array<int,2> temp1(1,3,FortranArray<2>());
    temp1=0;
    Array<int,2> temp2(1,3,FortranArray<2>());
    temp1=0;
    
    temp1(1,1)=seq1(1,i-2);
    temp1(1,2)=seq1(1,i-1);
    temp1(1,3)=seq1(1,i);
    
    temp2(1,1)=seq1(1,i-2);
    temp2(1,2)=seq1(1,i-1);
    temp2(1,3)=seq2(DAY,i);
    
    int code_old=0;
    int code_new=0;
    for (int j=1; j<=62; j++)
    {
        if (temp1(1,1)==codon(j,1)&&temp1(1,2)==codon(j,2)&&temp1(1,3)==codon(j,3))
        {
            code_old=codon(j,4);
        }
        if (temp2(1,1)==codon(j,1)&&temp2(1,2)==codon(j,2)&&temp2(1,3)==codon(j,3))
        {
            code_new=codon(j,4);
        }
    }
    if (code_new==code_old)
    {
        return 1; //S
    }
    else
    {
        return 0;//NS
    }
    
}
#endif /* NS_judge_hpp */
