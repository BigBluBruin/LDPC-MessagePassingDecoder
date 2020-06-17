#include <vector>
#include <numeric>
#include <cmath>
#include <math.h>
#include <algorithm>
#include <iostream>

//--------------------------------------------------
double sgn(double  input);
double vari_node_operation(std::vector<double> input);
void LP(double &in, int l, int r);
//---------------------------------------------------
double boxplus(double argument_1, double argument_2);
double check_min_double(double & input_1, double & input_2);
double minstar(double & input_1, double & input_2);
//---------------------------------------------------
double check_node_operation(std::vector<double> input);
double check_node_operation_min_double(std::vector<double> & input);
double check_node_operation_minstar(std::vector<double> & input);
double check_node_operation_att_min_double(std::vector<double> & input, const double att);
