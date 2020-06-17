#pragma once
#include "Parity_Check_Matrix_Info.h"
#include "full_BP.h"
#include <iostream>
#include <string.h>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <random>

class ldpc_full_precision_decoder{

public:
std::string h_filename;
std::vector<double> parameters; //noise
Parity_Check_Matrix_Info h_ins;
int target_error;


public:
std::vector<int> total_frames;
std::vector<int> total_iterations;
int max_iter;
char* check_op;



public:
ldpc_full_precision_decoder(std::string H_filename, std::vector<double> Parameters,int Target_error, int Max_ter, char Check_op[]);
void noise_generator(std::vector<double> & cwds,double parameter);
void fill_in();
bool iscwds(std::vector<int> final_bits);

//-----Decoder overloads , can be combined by one!-----
bool decoder(std::vector<double> cwds,int &iteration, std::vector<int> & final_bits, double sigma2);
bool decoder(std::vector<double> cwds,int &iteration, std::vector<int> & final_bits,std::vector<int> & initial_bits ,  double sigma2);
void decoder_track(std::vector<double> cwds,int &iteration, std::vector<int> & final_bits,std::vector<int> & initial_bits ,  double sigma2, int ind);


//-----main simulation overloads-----
void main_simulation(const char ind[]); // all zero codeword 
void main_simulation(const char ind[], std::string G_filename); // randomly generate codeword

//----- random codeword generation -----
void generate_whole_bits(std::vector<int> & inputbit);
void generate_codeword(std::vector<int> &inputbit, std::vector<double> & codewords);
};