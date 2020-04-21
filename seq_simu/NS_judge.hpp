//
//  NS_judge.hpp
//  seq_simu
//
//  Created by Taolee on 3/29/20.
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
    Array<int,2> temp1(3,1,FortranArray<2>());
    temp1=0;
    Array<int,2> temp2(3,1,FortranArray<2>());
    temp1=0;
    
    temp1(1,1)=seq1(i-2,1);
    temp1(2,1)=seq1(i-1,1);
    temp1(3,1)=seq1(i,1);
    
    temp2(1,1)=seq1(i-2,1);
    temp2(2,1)=seq1(i-1,1);
    temp2(3,1)=seq2(i,DAY);
    
    int code_old=1;
    int code_new=0;
    for (int j=1; j<=62; j++)
    {
        if (temp1(1,1)==codon(1,j)&&temp1(2,1)==codon(2,j)&&temp1(3,1)==codon(3,j))
        {
            code_old=codon(4,j);
        }
        if (temp2(1,1)==codon(1,j)&&temp2(2,1)==codon(2,j)&&temp2(3,1)==codon(3,j))
        {
            code_new=codon(4,j);
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
