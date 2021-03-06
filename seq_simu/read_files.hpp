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
void read_file(Array<int, 2> &seq_orig_sar,Array<int, 2> &seq_orig_tg,Array<int, 2> &codon,Array<double, 2> &sitefreq,Array<double, 2> &sitefreq_out1, Array<double, 2> &sitefreq_out2, int seq_length, string evolving_sequence, string reference_sequence, string codon_file, string substitution_file, string substitution_file_outgroup1, string substitution_file_outgroup2)
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
        seq_orig_sar(i+1,1)=array[i];
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
        seq_orig_tg(i+1,1)=array1[i];
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
            codon(j+1,i+1)=syno[i][j];
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
            sitefreq(j+1,i+1)=site_freq[i][j];
        }
    }
    infile21.close();
    
    
    double site_freq_out1[12][4]={0};
       ifstream infile22;
       infile22.open(substitution_file_outgroup1);
       double* ptr22 = &site_freq_out1[0][0];
       while(!infile22.eof())
       {
           infile22>>*ptr22;
           ptr22++;
       }
       for (int i=0; i<12;i++)
       {
           for (int j=0;j<4;j++)
           {
               sitefreq_out1(j+1,i+1)=site_freq_out1[i][j];
           }
       }
       infile22.close();
    
    
     double site_freq_out2[12][4]={0};
          ifstream infile23;
          infile23.open(substitution_file_outgroup2);
          double* ptr23 = &site_freq_out2[0][0];
          while(!infile23.eof())
          {
              infile23>>*ptr23;
              ptr23++;
          }
          for (int i=0; i<12;i++)
          {
              for (int j=0;j<4;j++)
              {
                  sitefreq_out2(j+1,i+1)=site_freq_out2[i][j];
              }
          }
          infile23.close();
}

#endif /* read_files_hpp */
