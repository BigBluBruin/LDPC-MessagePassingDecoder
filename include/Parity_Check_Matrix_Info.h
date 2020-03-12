#pragma once
#include <vector>
#include <fstream>
#include <string>
#include <iostream>

//This class defines the information of an IB decoder
class Parity_Check_Matrix_Info
{
    public:
    std::string filename;


    public:
    int vari_num;
    int check_num;
    int edge_num;
    int vari_max;
    int check_max;
    double rate;
    std::vector<std::vector <int>>  edge_v;
    std::vector<std::vector <int>>  edge_c;
    std::vector <int>  vari_degreetable;
    std::vector <int> check_degreetable;
    std::vector<std::vector <int>>  edge_relation;
    std::vector<std::vector<unsigned>>  G_info;

    public:
    Parity_Check_Matrix_Info(std::string Filename="None");
    bool Read_Parity_Check_Matrix();
    bool Read_G_Matrix(std::string fiilename);
};
