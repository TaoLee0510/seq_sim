# seq_sim

Seq_sim is an in-silico evolutionary process to identify the confidence of a phylogeny tree based on a specific out-group or out-groups. To compile the codes, the blitz++ and gsl are required. The gsl may NOT compatible with the Windows system, please switch to macOS or a Linux based system.

To compile the codes, the blitz++ and gsl are required. The gsl may NOT compatible with the Windows system, please switch to macOS or a Linux based system.

Author: Li Tao

Need help? Contact me via e-mail: taolee.lit@gmail.com

Install:
For macOS user:
Import the project to Xcode to use it directly.

For Linux user:(clang is highly recommanded)

clang++ main.cpp -lgsl -lgslcblas -v -o seq_sim

OR

g++ main.cpp -lgsl -lgslcblas -v -o seq_sim

Usage:

./seq_sim --evolving_sequence=path --reference_sequence=path --codon_file=path --substitution_file=path --substitution_file_outgroup1=path --substitution_file_outgroup2=path --seq_length=int --mutation_rate=double --Evolution_days=int --Pre_evolution_days=int --time_interval=int --simulation_times=int --Divergence=double --Divergency_sampling_times=int --output_evo_seq_file=int --KaKS=double --n_selection=double

example:

./seq_simu --evolving_sequence=./Sars_cov_2_orig.txt --reference_sequence=./TG13_orig.txt --codon_file=./ns_and_s.txt --substitution_file=./site_freq.txt --substitution_file_outgroup1=./site_freq_outgroup1.txt --substitution_file_outgroup2=./site_freq_outgroup2.txt --seq_length=29274 --mutation_rate=0.2542 --mutation_rate_outgroup1=0.2537 --mutation_rate_outgroup2=0.2534 --Evolution_days=30000 --Pre_evolution_days=27000 --time_interval=100 --simulation_times=34 --Divergence=0.005 --Divergency_sampling_times=1000 --output_evo_seq_file=0 --KaKS=0.05 --n_selection=0
  
  Input:
  
  evolving_sequence: The target sequence to be evolving. Any sequence is suitable. 
  
  reference_sequence: A reference sequence to calculate the similarity at the whole sequence level.
  
  codon_file: A file with 4 columns include the codons and the amino acid coded by a Condon. The first 3 columns represent the codon and the last column represents the amino acid.
  
  substitution_file: A file with 4 columns include the base substitution rate.  The first columns represent the original base, the second is the substitution base, the third the substitution rate. The last column represents the cumulative probability (rank the substitution rate from small to large). 

  seq_length: The sequence length of an aligned evolving_sequence.
  
  mutation_rate: Mutation rate. unit: sites/day.
  
  Evolution_days: Days of the evolving_sequence evolved.
  
  Pre_evolution_days: Additional evolution times of sequences to generate a "outgroup" with a further distance.
  
  time_interval: This parameter sets the time interval to record the evolved sequence. Meanwhile, at the recording time, the program calculated confidence.
  
  simulation_times: The times of the simulations you wanna run.
  
  Divergence: This parameter controls the polymorphic level of two sequences which are used to calculate the phylogenic confidence.
  
  Divergency_sampling_times: The resampling time of the polymorphic sequences.
  
  output_evo_seq_file: A int type parameter to control the output files whether include the evolved sequences. 0: DO NOT output. 1: output.
  
  KaKS: The Value of Ka/Ks.
  
  n_selection: A Ka/Ks value to control the evolutionary process during t3.
    
  Note:
  
  1, In order to make the program run easily, all Base and amino acids names should convert into digital. Here we stipulate the A=1, T=2, G=3, C=4, others (including gaps)=0. the amino acids could convert into any int type digital as you like.
  
  2, The formate of evolving_sequence and reference_sequence is a txt formate file with space character is used as the column delimiter.
  
  3, The codon_file is a 62 x 4 (row x column) file (codon.txt). The substitution_file is a 12 x 4 (row x column) file (substitution.txt). Users can generate these two files to suit their objects. if the dimension is not the same as the examples, one easy way is to fill the related rows with 0. OR recreated the dimension settings on related .hpp files.
        
  Output:
  
  The program output .txt files:
  
  1, Parameters.txt: list the input parameters to run the program.
  
  
  2, similarity_ref_evo.txt: The divergency between ref and evo sequences. [colunms: filelables; over_all_divergency; synonymous_site_divergency; non-synonymous_site_divergency].
  
  
  3, Similarity_(n).txt: The divergency between N1 and N2 of nth simulation. [colunms: days(value x time_interval); over_all_divergency; synonymous_site_divergency; non-synonymous_site_divergency].
  
  
  4, Similarity_divergency_(n).txt: The divergency between N3 and N4 of nth simulation. [colunms: days(value x time_interval); over_all_divergency; synonymous_site_divergency; non-synonymous_site_divergency].
  
  
  5, results_single_out_(n).txt: Results of N1 as outgroup. [$1:over_all_correct_number; $2:over_all_error_number; $3: over_all_error_rate; $4: synonymous site_correct_number; $5: synonymous site_error_number; $6: synonymous site_error_rate; $7: non-synonymous_correct_number; $8: non-synonymous_error_number; $9: non-synonymous site_error_rate; $10: uncertainty_number; $11: uncertainty_rate; $12: synonymous uncertainty_number; $13: synonymous uncertainty_rate; $14: non-synonymous uncertainty_number; $15: non-synonymous uncertainty_rate].
  
  
  6, results_double_out_(n).txt: Results of N1 and N6 as outgroups. [$1:over_all_correct_number; $2:over_all_error_number; $3: over_all_error_rate; $4: synonymous site_correct_number; $5: synonymous site_error_number; $6: synonymous site_error_rate; $7: non-synonymous_correct_number; $8: non-synonymous_error_number; $9: non-synonymous site_error_rate; $10: uncertainty_number; $11: uncertainty_rate; $12: synonymous uncertainty_number; $13: synonymous uncertainty_rate; $14: non-synonymous uncertainty_number; $15: non-synonymous uncertainty_rate].
  
  
   7, mutation.txt: A mutation number record during time. [$1:Total mutation; $2-$5: inner check number, omit them; $6: mutaions ocurred during whole simulation time; $7: synonymous mutations ocurred during whole simulation time; $8: Non-synonymous mutations ocurred during whole simulation time; $9: Total mutations fixed; $10: Synonymous mutations fixed; $11: Non-synoymous mutations fixed]
   
  
  If the output_evo_seq_file=1, the following files would be outputed.
  
  7, Seq_evo_(n).txt: N2; each row represnts days(value x time_interval); n: nth simulation.
  
  8, Seq_evo1_(n).txt:N1; each row represnts days(value x time_interval); n: nth simulation.
  
  9, Seq_evo2_(n).txt:N6; each row represnts days(value x time_interval); n: nth simulation.
  
  
  
  
  
  
