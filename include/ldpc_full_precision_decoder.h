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

private:
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
bool decoder(std::vector<double> cwds,int &iteration, std::vector<int> & final_bits, double sigma2);
bool iscwds(std::vector<int> final_bits);
void main_simulation(const char ind[]);
};