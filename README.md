# seq_sim
in-silico evolutionary process to identify the confidence of phylogeny tree based on a specific out-group or out-groups. 
To complie the codes, the blitz++ and gsl are required. The gsl may NOT compatible with Windows system, please switch to MacOS or a Linux based system.

Need help? Contact me via e-mail: aeolia.zafiro@gmail.com

Usage:

./seq_sim --evolving_sequence=path --reference_sequence=path --codon_file=path --substitution_file=path --seq_length=int --mutation_rate=double --Evolution_days=int --Pre_evolution_days=int --time_interval=int --simulation_times=int --Divergence=double --Divergency_sampling_times=int --output_evo_seq_file=int --KaKS=double

example:

./seq_simu --evolving_sequence=./Sars_cov_2_orig.txt --reference_sequence=./TG13_orig.txt --codon_file=./codon.txt --substitution_file=./substitution.txt --seq_length=29274 --mutation_rate=0.2542 --Evolution_days=30000 --Pre_evolution_days=27000 --time_interval=100 --simulation_times=200 --Divergence=0.001 --Divergency_sampling_times=1000 --output_evo_seq_file=0 --KaKS=0.05
  
  Input:
  
  evolving_sequence: The target sequence to be evolving. Any sequence is suitable. 
  
  reference_sequence: A reference sequence to calculate the similarity at the whole sequence level.
  
  codon_file: A file with 4 columns include the codons and the amino acid coded by a Condon. The first 3 columns represent the codon and the last column represents the amino acid.
  
  substitution_file: A file with 4 columns include the base substitution rate.  The first columns represent the original base, the second is the substitution base, the third the substitution rate. The last column represents the cumulative probability (rank the substitution rate from small to large). 

  seq_length: The sequence length of an aligned evolving_sequence.
  
  mutation_rate: Mutation rate. unit: sites/day.
  
  Evolution_days: Days of the evolving_sequence evolved.
  
  Pre_evolution_days: Additional evolution times of a sequence to generate a "outgroup" with a further distance.
  
  time_interval: This parameter sets the time interval to record the evolved sequence. Meanwhile, at the recording time, the program calculated confidence.
  
  simulation_times: The times of the simulations you wanna run.
  
  Divergence: This parameter controls the polymorphic level of two sequences which are used to calculate the phylogenic confidence.
  
  Divergency_sampling_times: The resampling time of the polymorphic sequences.
  
  output_evo_seq_file: A int type parameter to control the output files whether include the evolved sequences. 0: DO NOT output. 1: output.
  
  KaKS: The Value of Ka/Ks.
  
    
  Note:
  
  1, In order to make the program run easily, all Base and amino acids names should convert into digital. Here we stipulate the A=1, T=2, G=3, C=4, others (including gaps)=0. the amino acids could convert any int type digital as you like.
  
  2, The evolving_sequence reference_sequence should be aligned before inputting to the program.
  
  3, The formate of evolving_sequence and reference_sequence is a txt formate file with space character is used as the column delimiter.
  
  4, The codon_file is a 62 x 4 (row x column) file (codon.txt). The substitution_file is a 12 x 4 (row x column) file (substitution.txt). Users can generate these two files to suit their objects. if the dimension is not the same as the examples, one easy way is to fill the related rows with 0. OR recreated the dimension settings on related .hpp files.
        
  Output:
  
  The program output .txt files:
  
  1, Parameters.txt: list the input parameters to run the program.
  
  2, similarity_ref_evo.txt: The divergency between ref and evo sequences. [colunms: filelables; over_all_divergency; synonymous_site_divergency; non-synonymous_site_divergency].
  
  3, Similarity_(n).txt: The divergency between N1 and N2 of nth simulation. [colunms: days(value x time_interval); over_all_divergency; synonymous_site_divergency; non-synonymous_site_divergency].
  
  4, Similarity_divergency_(n).txt: The divergency between N3 and N4 of nth simulation. [colunms: days(value x time_interval); over_all_divergency; synonymous_site_divergency; non-synonymous_site_divergency].
  
  5, results_single_out_(n).txt: Results of N1 as outgroup. [colunms: over_all_correct_number; over_all_error_number; over_all_error_rate; synonymous site_correct_number; synonymous site_error_number; synonymous site_error_rate; non-synonymous_correct_number; non-synonymous_error_number; non-synonymous site_error_rate; uncertainty_number; uncertainty_rate].
  
  6, results_double_out_(n).txt: Results of N1 and N6 as outgroups. [colunms: over_all_correct_number; over_all_error_number; over_all_error_rate; synonymous site_correct_number; synonymous site_error_number; synonymous site_error_rate; non-synonymous_correct_number; non-synonymous_error_number; non-synonymous site_error_rate; uncertainty_number; uncertainty_rate].
  
  If the output_evo_seq_file=1, the following files would be outputed.
  
  7, Seq_evo_(n).txt: N2; each row represnts days(value x time_interval); n: nth simulation.
  
  8, Seq_evo1_(n).txt:N1; each row represnts days(value x time_interval); n: nth simulation.
  
  9, Seq_evo2_(n).txt:N6; each row represnts days(value x time_interval); n: nth simulation.
  
  
  
  
  
  
