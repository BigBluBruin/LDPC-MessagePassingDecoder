#include "ldpc_full_precision_decoder.h"

ldpc_full_precision_decoder::ldpc_full_precision_decoder(std::string H_filename, std::vector<double> Parameters,int Target_error, int Max_iter, char Check_op[])
{
    h_filename=H_filename;
    parameters=Parameters;
    target_error=Target_error;
    max_iter=Max_iter;
    check_op=Check_op;
}

void ldpc_full_precision_decoder::noise_generator(std::vector<double> &cwds, double parameter)
{
    std::random_device rd;
    std::mt19937 generator(rd());
    std::normal_distribution<double> distribution(0, sqrt(parameter));
    for (unsigned inc = 0; inc < cwds.size(); inc++)
    {
        cwds[inc] += distribution(generator);
    }
}

void ldpc_full_precision_decoder::fill_in()
{
    h_ins.filename = h_filename;
    if (h_ins.Read_Parity_Check_Matrix())
    {
        std::cout << "Parity Matrix Success"
                  << std::endl;
    }
    else 
    {
        std::cout << "Parity Check Matrix Fails"
                  << std::endl;
    }
}

bool ldpc_full_precision_decoder::decoder(std::vector<double> cwds, int &iteration, std::vector<int> &final_bits, double sigma2)
{
    //Definition area
    std::vector<double> msg_c2v(h_ins.edge_num, -1);
    std::vector<double> msg_v2c(h_ins.edge_num, -1);
    std::vector<double> rx(h_ins.vari_num, -1);
    double final_codewords;
    std::vector<double> message;


    for (int ii = 0; ii < h_ins.vari_num; ii++)
    {
        rx[ii] = (2/sigma2)*cwds[ii];
    }

    //Ini v2c messages
    for (int ii = 0; ii < h_ins.vari_num; ii++)
    {
        for (int jj = 0; jj < h_ins.vari_degreetable[ii]; jj++)
        {
            msg_v2c[h_ins.edge_v[ii][jj]] = rx[ii];
        }
    }

    //Start Iteration
    for (int cur_iter = 0; cur_iter < max_iter; cur_iter++)
    {
        //c2v update
        for (int ii = 0; ii < h_ins.check_num; ii++)
        {
            int cur_dc = h_ins.check_degreetable[ii];
            for (int jj = 0; jj < cur_dc; jj++)
            {
                message.clear();
                for (int kk = 0; kk < cur_dc; kk++)
                {
                    if (kk != jj)
                    {
                        message.push_back(msg_v2c[h_ins.edge_c[ii][kk]]);
                    }
                }
                if(strcmp(check_op,"boxplus")==0)
                {
                    msg_c2v[h_ins.edge_c[ii][jj]]=check_node_operation(message);
                }
                else if (strcmp(check_op,"minsum")==0)
                {
                    msg_c2v[h_ins.edge_c[ii][jj]]=check_node_operation_min_double(message);
                }
                else if (strcmp(check_op,"att_minsum")==0)
                {
                    msg_c2v[h_ins.edge_c[ii][jj]]=check_node_operation_att_min_double(message,0.7);
                }
                else if (strcmp(check_op,"minstar")==0)
                {
                    msg_c2v[h_ins.edge_c[ii][jj]]=check_node_operation_minstar(message);
                }
                else
                {
                    std::cout<<"Info: no such operation "<< check_op<<".... please check again"<<std::endl;
                }
            }
        }

        //v2c update
        for (int ii = 0; ii < h_ins.vari_num; ii++)
        {
            int cur_dv = h_ins.vari_degreetable[ii]; //find current dv
            for (int jj = 0; jj < cur_dv; jj++)
            {
                message.clear();
                message.push_back(rx[ii]);
                //std::cout<<"have been 2.5"<<std::endl;
                for (int kk = 0; kk < cur_dv; kk++)
                {
                    //collect data
                    if (kk != jj)
                    {
                        message.push_back(msg_c2v[h_ins.edge_v[ii][kk]]);
                    }
                }
                msg_v2c[h_ins.edge_v[ii][jj]] = vari_node_operation(message);
            }
        }

        //final decision
        for (int ii = 0; ii < h_ins.vari_num; ii++)
        {

            int cur_dv = h_ins.vari_degreetable[ii];
            message.clear();
            message.push_back(rx[ii]);
            for (int jj = 0; jj < cur_dv; jj++)
            {
                message.push_back(msg_c2v[h_ins.edge_v[ii][jj]]);
            }
            final_codewords = vari_node_operation(message);
            if (final_codewords > 0)
            {
                final_bits[ii] = 0;
            }
            else
            {
                final_bits[ii] = 1;
            }
        }

        //check sum
        if (iscwds(final_bits))
        {
            iteration = iteration + cur_iter + 1;
            int sum = 0;
            for (int ii = 0; ii < h_ins.vari_num; ii++)
            {
                sum = sum + final_bits[ii];
            }
            if (sum == 0)
            {
                return true;
            }
            else
            {
                return false;
            }
        }
    }
    iteration = iteration + max_iter;
    return false;
}

bool ldpc_full_precision_decoder::iscwds(std::vector<int> final_bits)
{
    std::vector<int> edge_bits(h_ins.edge_num, -1);
    for (int ii = 0; ii < h_ins.vari_num; ii++)
    {
        for (int jj = 0; jj < h_ins.vari_degreetable[ii]; jj++)
        {
            edge_bits[h_ins.edge_v[ii][jj]] = final_bits[ii];
        }
    }

    for (int ii = 0; ii < h_ins.check_num; ii++)
    {
        int checks = 0;
        for (int jj = 0; jj < h_ins.check_degreetable[ii]; jj++)
        {
            checks = (checks + edge_bits[h_ins.edge_c[ii][jj]]) % 2;
        }
        if (checks > 0)
        {
            return false;
        }
    }

    return true;
}

void ldpc_full_precision_decoder::main_simulation(const char ind[])
{
    int total_num = 0;
    int error_num = 0;
    int iteration;
    double cur_para;
    std::vector<double> codewords;
    std::vector<int> final_bits;
    for (unsigned ii = 0; ii < parameters.size(); ii++)
    {
        iteration = 0;
        cur_para = pow(10, (-0.1 * parameters[ii]) / (2.0 * h_ins.rate));
        std::string result_file = "Result_n_" + std::to_string(h_ins.vari_num) + "_k_" + std::to_string(h_ins.check_num) + "_check_" + check_op + "_Para_" + std::to_string(parameters[ii]) + "_ind_" + ind + ".txt";
        std::ofstream result_bar;
        do
        {
            total_num = total_num + 1;
            codewords.assign(h_ins.vari_num, 1);
            final_bits.assign(h_ins.vari_num, -1);
            //add noise
            noise_generator(codewords, cur_para);
            bool result = decoder(codewords, iteration, final_bits,cur_para);
            if (!result)
            {
                error_num = error_num + 1;
                std::cout << "Para: " << parameters[ii] << " error: " << error_num << "fer: " << (double)error_num / (double)total_num << std::endl;
            }
            if (total_num % 100 == 1)
            {
                result_bar.open(result_file);
                result_bar << parameters[ii] << "  " << total_num << "  " << error_num << "  " << (double)iteration / (double)total_num << std::endl;
                result_bar.close();
            }

        } while (error_num < target_error);
        result_bar.open(result_file);
        result_bar << parameters[ii] << "  " << total_num << "  " << error_num << "  " << (double)iteration / (double)total_num << std::endl;
        result_bar.close();
    }
}