//
//  read_files.hpp
//  seq_simu
//
//  Created by Taolee on 3/30/20.
//  Copyright © 2020 Taolee. All rights reserved.
//

#ifndef read_files_hpp
#define read_files_hpp

#include <stdio.h>
#include <iostream>
#include <string>
using namespace std;
using namespace blitz;
void read_file(Array<int, 2> &seq_orig_sar,Array<int, 2> &seq_orig_tg,Array<int, 2> &codon,Array<double, 2> &sitefreq,int seq_length, string evolving_sequence, string reference_sequence, string codon_file, string substitution_file)
{
    int array[29274]={0};
    ifstream infile;
    infile.open(evolving_sequence);
    int* ptr = &array[0];
    while(!infile.eof())
    {
        infile>>*ptr;
        ptr++;
    }
    for (int i=0; i<seq_length;i++)
    {
        seq_orig_sar(1,i+1)=array[i];
    }
    infile.close();

    int array1[29274]={0};
    ifstream infile1;
    infile1.open(reference_sequence);
    int* ptr1 = &array1[0];
    while(!infile1.eof())
    {
        infile1>>*ptr1;
        ptr1++;
    }
    for (int i=0; i<seq_length;i++)
    {
        seq_orig_tg(1,i+1)=array1[i];
    }
    infile1.close();


    int syno[62][4]={0};
    ifstream infile11;
    infile11.open(codon_file);
    int* ptr11 = &syno[0][0];
    while(!infile11.eof())
    {
        infile11>>*ptr11;
        ptr11++;
    }
    for (int i=0; i<62;i++)
    {
        for (int j=0;j<4;j++)
        {
            codon(i+1,j+1)=syno[i][j];
        }
    }
    infile11.close();
    
    double site_freq[12][4]={0};
    ifstream infile21;
    infile21.open(substitution_file);
    double* ptr21 = &site_freq[0][0];
    while(!infile21.eof())
    {
        infile21>>*ptr21;
        ptr21++;
    }
    for (int i=0; i<12;i++)
    {
        for (int j=0;j<4;j++)
        {
            sitefreq(i+1,j+1)=site_freq[i][j];
        }
    }
    infile21.close();
}

#endif /* read_files_hpp */
