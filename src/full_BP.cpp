#include "full_BP.h"

double boxplus(double argument_1, double argument_2)
{
    double result = 2 * atanh(tanh(argument_1 / 2) * tanh(argument_2 / 2));
    if (result==INFINITY)
    {
        result = 3.0*pow(10,7);
    }
    else if (result==-INFINITY)
    {
        result = - 3.0*pow(10,7);
        //std::cout<<"here"<<std::endl;
    }
    return result;
}

double check_node_operation(std::vector<double> input)
{
    double output, fir, sec;
    fir = input[0];
    for (unsigned index = 0; index < input.size() - 1; index++)
    {
        sec = input[index + 1];
        fir = boxplus(fir, sec);
    }
    output = fir;
    return output;
}

double vari_node_operation(std::vector<double> input)
{
    double sum;
    sum = std::accumulate(input.begin(), input.end(), 0.0);
    return sum;
}



double sgn(double input)
{
    if (input >= 0)
    {
        return 1;
    }
    else
    {
        return -1;
    }
}
double check_min_double(double &input_1, double &input_2)
{
    return sgn(input_1) * sgn(input_2) * std::min(abs(input_2), abs(input_1));
}

double check_node_operation_min_double(std::vector<double> &input)
{
    double fir = input[0];
    double sec;
    for (unsigned ii = 0; ii < input.size() - 1; ii++)
    {
        sec = input[ii + 1];
        fir = check_min_double(fir, sec);
    }
    return fir;
}

double check_node_operation_att_min_double(std::vector<double> &input, const double att)
{
    double fir = input[0];
    double sec;
    for (unsigned ii = 0; ii < input.size() - 1; ii++)
    {
        sec = input[ii + 1];
        fir = check_min_double(fir, sec);
    }
    return att * fir;
}

double minstar(double &input_1, double &input_2)
{
    double part1= -sgn(input_1)*sgn(input_2)*std::min(abs(input_1), abs(input_2));
    double part2= log(1+exp(-abs(input_2-input_1)));
    double part3= log(1+exp(-abs(input_1+input_2)));
    return sgn(input_1)*sgn(input_2)*-1*(part1+part2-part3);
}

double check_node_operation_minstar(std::vector<double> & input)
{
    double fir = input[0];
    double sec;
    for (unsigned ii = 0; ii < input.size() - 1; ii++)
    {
        sec = input[ii + 1];
        fir = minstar(fir, sec);
    }
    return fir;
}

void LP(double &in, int l, int r)
{

    double rpow, lpow, minpow, maxpow;
    rpow = pow(2, r);
    lpow = pow(2, (l - 1));
    minpow = (-1) * lpow;
    maxpow = lpow - (double)1.0 / (double)rpow;
    if (in > maxpow)
    {
        (in) = maxpow;
    }
    else if (in < minpow)
    {
        (in) = minpow;
    }
    else
    {
        (in) = (double)floor((in)*rpow + 0.5) / (double)rpow;
    }
}