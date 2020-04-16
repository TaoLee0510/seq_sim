//
//  similarity_tg.hpp
//  seq_simu
//
//  Created by Taolee on 3/29/20.
//  Copyright © 2020 Taolee. All rights reserved.
//

#ifndef similarity_tg_hpp
#define similarity_tg_hpp

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
void similarity_tg(Array<int, 2> codon, Array<int, 2> seq_orig_sar,Array<int, 2> seq_orig_tg,Array<double, 2> &Similarity_tg,int seq_length,int times)
{
    int NS_similarity=0;
    int S_similarity=0;
    for (int i=1;i<=seq_length;i++)
    {
        if (seq_orig_sar(i,1)!=seq_orig_tg(i,1))
        {
            if (i%3==0)
            {
                int s=NS_judge(i,seq_orig_sar,seq_orig_tg,codon,1);
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
    }
    Similarity_tg(times,1)=1;
    double NS_site_number=seq_length*3.5/4.5;
    double S_site_number=seq_length/4.5;
    double sum5=S_similarity+NS_similarity;
    Similarity_tg(times,2)=sum5/seq_length;
    Similarity_tg(times,3)=S_similarity/S_site_number;
    Similarity_tg(times,4)=NS_similarity/NS_site_number;
}
#endif /* similarity_tg_hpp */
