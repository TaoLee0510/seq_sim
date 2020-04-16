//
//  main.cpp
//  seq_simu Ver. 5.0
//
//  Created by Taolee on 3/29/20.
//  Copyright © 2020 Taolee. All rights reserved.
//

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
#include <string>
#include "Seq_simu.hpp"

using namespace std;
#pragma pack(8)
static const char *short_options = "E:R:C:S:l:u:D:d:i:T:v:s:F";
static const struct option long_options[] = {
    {"evolving_sequence", required_argument, NULL, 'E'},
    {"reference_sequence", required_argument, NULL, 'R'},
    {"codon_file", required_argument, NULL, 'C'},
    {"substitution_file", required_argument, NULL, 'S'},
    {"seq_length", required_argument, NULL, 'l'},
    {"mutation_rate", required_argument, NULL, 'u'},
    {"Evolution_days", required_argument, NULL, 'D'},
    {"Pre_evolution_days", required_argument, NULL, 'd'},
    {"time_interval", required_argument, NULL, 'i'},
    {"simulation_times", required_argument, NULL, 'T'},
    {"Divergency_mutations", required_argument, NULL, 'v'},
    {"Divergency_sampling_times", required_argument, NULL, 's'},
    {"output_evo_seq_file", required_argument, NULL, 'F'},
    {NULL, 0, NULL, 0}
};
int main(int argc, char *argv[]) {
    string evolving_sequence;
    string reference_sequence;
    string codon_file;
    string substitution_file;
    int seq_length = 0;
    double mutation_rate = 0.0;
    int Evolution_days = 0;
    int Pre_evolution_days = 0;
    int time_interval = 0;
    int simulation_times = 0;
    int Divergency_mutations = 0;
    int Divergency_sampling_times =0;
    int output_evo_seq_file=0;//0 no; 1 yes;
    int opt = 0;
    while( (opt = getopt_long(argc, argv, short_options, long_options, NULL)) != -1){
        switch (opt){
            case '?':
                fprintf(stdout, "Usage: %s --evolving_sequence=<path> --reference_sequence=<path> --codon_file=<path> --substitution_file=<path> --seq_length=<int> --mutation_rate=<double> --Evolution_days=<int> --Pre_evolution_days=<int> --time_interval=<int> --simulation_times=<int> --Divergency_mutations=<int> --Divergency_sampling_times=<int>", argv[0]);
                return 0;
            case 'E':
                evolving_sequence = optarg;
                break;
            case 'R':
                reference_sequence = optarg;
                break;
            case 'C':
                codon_file = optarg;
                break;
            case 'S':
                substitution_file = optarg;
                break;
            case 'l':
                seq_length = atoi(optarg);
                break;
            case 'u':
                mutation_rate = atof(optarg);
                break;
            case 'D':
                Evolution_days = atoi(optarg);
                break;
            case 'd':
                Pre_evolution_days = atoi(optarg);
                break;
            case 'i':
                time_interval = atoi(optarg);
                break;
            case 'T':
                simulation_times = atoi(optarg);
                break;
            case 'v':
                Divergency_mutations = atoi(optarg);
                break;
            case 's':
                Divergency_sampling_times = atoi(optarg);
                break;
            case 'F':
                output_evo_seq_file = atoi(optarg);
                break;
        }
    }
    
    //////////////////////////////////////////////////// paramters initiation/////////////////////////////////////////////
    double probability_of_mutation_sar=mutation_rate;
    int DDAAYY=Evolution_days;
    int DAY_pa=Pre_evolution_days;
    int interval=time_interval;
    int times=simulation_times;
    int D=DDAAYY/interval;
    int random_mutation_number=Divergency_mutations;
    //////////////////////////////////////////////////////////////////////////output parameters/////////////////////////////////////////////////////
    char filedir [100] = {'\0'};
    sprintf(filedir, "./Parameters.txt");
    FILE * fid1;
    fid1=fopen (filedir,"w+");
    fprintf(fid1, "%s %s %s\n" ,"Evolving Sequence", "=", evolving_sequence.c_str());
    fprintf(fid1, "%s %s %s\n" ,"Reference Sequence", "=", reference_sequence.c_str());
    fprintf(fid1, "%s %s %s\n" ,"Codon file", "=", codon_file.c_str());
    fprintf(fid1, "%s %s %s\n" ,"Substitution file", "=", substitution_file.c_str());
    fprintf(fid1, "%s %s %d\n" ,"Sequence length", "=", seq_length);
    fprintf(fid1, "%s %s %lf\n" ,"Mutation Rate", "=", mutation_rate);
    fprintf(fid1, "%s %s %d\n" ,"Evolution days", "=", Evolution_days);
    fprintf(fid1, "%s %s %d\n" ,"Preevolution days", "=", Pre_evolution_days);
    fprintf(fid1, "%s %s %d\n" ,"Time interval", "=", time_interval);
    fprintf(fid1, "%s %s %d\n" ,"simulation Times", "=", simulation_times);
    fprintf(fid1, "%s %s %d\n" ,"Divergency mutations", "=", Divergency_mutations);
    fprintf(fid1, "%s %s %d\n" ,"Divergency sampling times", "=", Divergency_sampling_times);
    fprintf(fid1, "%s %s %d\n" ,"Output evo_seq file", "=", output_evo_seq_file);
    fclose(fid1);
    ////////////////////////////////////////////////////////////////////////// run seq_simu /////////////////////////////////////////////////////
    Seq_simu(seq_length, probability_of_mutation_sar,DDAAYY,DAY_pa, interval, times, D, random_mutation_number,Divergency_sampling_times,evolving_sequence
             ,reference_sequence,codon_file,substitution_file,output_evo_seq_file);
    return 0;
}
