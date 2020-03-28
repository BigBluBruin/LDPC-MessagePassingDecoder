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
    std::vector<int> ini_bits(h_ins.vari_num,0);
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
            bool result = decoder(codewords, iteration, final_bits,ini_bits,cur_para);
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

void ldpc_full_precision_decoder::main_simulation(const char ind[], std::string G_filename)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::binomial_distribution<> d(1, 0.5);
    h_ins.Read_G_Matrix(G_filename);
    int total_num = 0;
    int error_num = 0;
    int iteration;
    int iteration_2=0;
    double cur_para;
    std::vector<double> codewords;
    std::vector<int> final_bits;
    std::vector<int> initial_bits;
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
            initial_bits.assign(h_ins.vari_num, 0);
            for (int ii = 0; ii < (h_ins.vari_num - h_ins.check_num); ii++)
            {
                initial_bits[ii] = d(gen);
            }
            //add noise
            generate_whole_bits(initial_bits);
            generate_codeword(initial_bits, codewords);
            noise_generator(codewords, cur_para);
            bool result = decoder(codewords, iteration, final_bits,initial_bits, cur_para);
            if (!result)
            {
                error_num = error_num + 1;
                int cur_ind = target_error*std::stoi(ind)+(error_num-1);
                std::string wcwd_filename="wrong_cwd_"+std::to_string(cur_ind)+".txt";
                std::string inib_filename="initial_bits_"+std::to_string(cur_ind)+".txt";
                std::ofstream wcwd_file(wcwd_filename);
                std::ofstream inib_file(inib_filename);
                for (unsigned ii =0; ii<codewords.size();ii++)
                {
                    wcwd_file<<codewords[ii]<<"  ";
                    inib_file<<initial_bits[ii]<<"  ";
                }
                wcwd_file.close();
                decoder_track(codewords,iteration_2,final_bits,initial_bits,cur_para,cur_ind);
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

void ldpc_full_precision_decoder::generate_whole_bits(std::vector<int> &inputbit)
{
    for (unsigned ii = 0; ii < h_ins.G_info[0].size(); ii++)
    {
        inputbit[h_ins.G_info[1][ii]] += inputbit[h_ins.G_info[0][ii]];
    }
    for (unsigned ii = 0; ii < inputbit.size(); ii++)
    {
        inputbit[ii] = inputbit[ii] % 2;
    }
}

void ldpc_full_precision_decoder::generate_codeword(std::vector<int> &inputbit, std::vector<double> &codewords)
{
    for (unsigned ii = 0; ii < inputbit.size(); ii++)
    {
        codewords[ii] = (inputbit[ii] == 1) ? -1.0 : 1.0;
    }
}

bool ldpc_full_precision_decoder::decoder(std::vector<double> cwds,int &iteration, std::vector<int> & final_bits,std::vector<int> & initial_bits ,  double sigma2)
{
    //Definition area
    std::vector<double> msg_c2v(h_ins.edge_num, -1);
    std::vector<double> msg_v2c(h_ins.edge_num, -1);
    std::vector<double> rx(h_ins.vari_num, -1);
    std::vector<double> forward, backward;
    double final_codewords;
    std::vector<double> message;
    int cur_dc,cur_dv;


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

    // //Start Iteration
    for (int cur_iter = 0; cur_iter < max_iter; cur_iter++)
    {


         //c2v update forward backword encoding
        for (int ii = 0; ii < h_ins.check_num; ii++)
        {
            cur_dc = h_ins.check_degreetable[ii];
            forward.assign(cur_dc, -1);
            backward.assign(cur_dc, -1);
            forward[0] = msg_v2c[h_ins.edge_c[ii][0]];
            for (int jj = 1; jj < cur_dc - 1; jj++)
            {
                if (strcmp(check_op, "boxplus") == 0)
                {
                    forward[jj] = boxplus(forward[jj - 1], msg_v2c[h_ins.edge_c[ii][jj]]);
                }
                else if (strcmp(check_op, "minsum") == 0)
                {
                    forward[jj] = check_min_double(forward[jj - 1], msg_v2c[h_ins.edge_c[ii][jj]]);
                }
                else
                {
                    std::cout << "Info: no such operation " << check_op << ".... please check again" << std::endl;
                }
            }
            msg_c2v[h_ins.edge_c[ii][cur_dc - 1]] = forward[cur_dc - 2];
            backward[cur_dc - 1] = msg_v2c[h_ins.edge_c[ii][cur_dc - 1]];
            for (int jj = cur_dc - 2; jj > 0; jj--)
            {
                if (strcmp(check_op, "boxplus") == 0)
                {
                    backward[jj] = boxplus(backward[jj + 1], msg_v2c[h_ins.edge_c[ii][jj]]);
                    msg_c2v[h_ins.edge_c[ii][jj]] = boxplus(backward[jj + 1], forward[jj - 1]);
                }
                else if (strcmp(check_op, "minsum") == 0)
                {
                    backward[jj] = check_min_double(backward[jj + 1], msg_v2c[h_ins.edge_c[ii][jj]]);
                    msg_c2v[h_ins.edge_c[ii][jj]] = check_min_double(backward[jj + 1], forward[jj - 1]);
                }
                else
                {
                    std::cout << "Info: no such operation " << check_op << ".... please check again" << std::endl;
                }
            }
            msg_c2v[h_ins.edge_c[ii][0]] = backward[1];
        }

        //v2c update
        for (int ii = 0; ii < h_ins.vari_num; ii++)
        {

            cur_dv = h_ins.vari_degreetable[ii];
            message.clear();
            message.push_back(rx[ii]);
            for (int jj = 0; jj < cur_dv; jj++)
            {
                message.push_back(msg_c2v[h_ins.edge_v[ii][jj]]);
            }
            final_codewords = vari_node_operation(message);
            for (int jj = 0; jj < cur_dv; jj++)
            {
                msg_v2c[h_ins.edge_v[ii][jj]] = final_codewords-msg_c2v[h_ins.edge_v[ii][jj]];
            }
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
            for (int ii = 0; ii < h_ins.vari_num; ii++)
            {
                if (final_bits[ii] != initial_bits[ii])
                {
                    return false;
                }
            }
            return true;
        }
    }
    iteration = iteration + max_iter;
    return false;
}



void ldpc_full_precision_decoder::decoder_track(std::vector<double> cwds, int &iteration, std::vector<int> &final_bits, std::vector<int> &initial_bits, double sigma2, int ind)
{
    //Definition area
    std::vector<double> msg_c2v(h_ins.edge_num, -1);
    std::vector<double> msg_v2c(h_ins.edge_num, -1);
    std::vector<double> final_codewords(h_ins.vari_num,-1);
    std::vector<double> rx(h_ins.vari_num, -1);
    std::string ts_filename = "trapping_set_boxplus_"+std::to_string(ind)+".txt";
    std::string llr_filename = "wrong_llr_boxplus_"+std::to_string(ind)+".txt";
    std::ofstream myfile(ts_filename);
    std::ofstream my_llr_file(llr_filename);
    std::vector<int> cur_ts;
    std::vector<double> cur_llr;

    std::vector<double> message;

    for (int ii = 0; ii < h_ins.vari_num; ii++)
    {
        rx[ii] = (2 / sigma2) * cwds[ii];
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
                if (strcmp(check_op, "boxplus") == 0)
                {
                    msg_c2v[h_ins.edge_c[ii][jj]] = check_node_operation(message);
                }
                else if (strcmp(check_op, "minsum") == 0)
                {
                    msg_c2v[h_ins.edge_c[ii][jj]] = check_node_operation_min_double(message);
                }
                else if (strcmp(check_op, "att_minsum") == 0)
                {
                    msg_c2v[h_ins.edge_c[ii][jj]] = check_node_operation_att_min_double(message, 0.7);
                }
                else if (strcmp(check_op, "minstar") == 0)
                {
                    msg_c2v[h_ins.edge_c[ii][jj]] = check_node_operation_minstar(message);
                }
                else
                {
                    std::cout << "Info: no such operation " << check_op << ".... please check again" << std::endl;
                }
                if (msg_c2v[h_ins.edge_c[ii][jj]] == -INFINITY || msg_c2v[h_ins.edge_c[ii][jj]] == INFINITY)
                {
                    for (auto aa : message)
                        std::cout << aa << " ";
                    std::cout << std::endl;
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
            final_codewords[ii] = vari_node_operation(message);

            if (final_codewords[ii] > 0)
            {
                final_bits[ii] = 0;
            }
            else
            {
                final_bits[ii] = 1;
            }
        }
        cur_ts.clear();
        cur_llr.clear();
        //check sum
        for (int ii = 0; ii < h_ins.vari_num; ii++)
        {
            if (final_bits[ii] != initial_bits[ii])
            {
                cur_ts.push_back(ii);
                cur_llr.push_back(final_codewords[ii]);
            }
        }

        //write out 
        for(unsigned ii=0;ii<cur_ts.size();ii++)
        {
            myfile<<cur_ts[ii]<<"  ";
            my_llr_file<<final_codewords[ii]<<" ";
        }
        myfile<<std::endl;
        my_llr_file<<std::endl;
        if(iscwds(final_bits))
            std::cout<<ind<<"is undetected cwds"<<std::endl;
    }
    iteration = iteration + max_iter;
    myfile.close();
    my_llr_file.close();
}



//-----trash code-----
/*cur_dv = h_ins.vari_degreetable[5];
for (int ii = 0; ii < cur_dc; ii++)
{
    std::cout << msg_c2v[h_ins.edge_v[5][ii]] << "  ";
}
std::cout << std::endl;

for (int ii = 0; ii < h_ins.vari_num; ii++)
{
    cur_dv = h_ins.vari_degreetable[ii]; //find current dv
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
}*/

/*cur_dv = h_ins.vari_degreetable[5];
for (int ii = 0; ii < cur_dc; ii++)
{
    std::cout << msg_c2v[h_ins.edge_v[5][ii]] << "  ";
}
std::cout << std::endl;
return false;

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

        cur_dc = h_ins.check_degreetable[0];
        for (int ii = 0; ii < cur_dc; ii++)
        {
            std::cout << msg_c2v[h_ins.edge_c[0][ii]] << "  ";
        }
        std::cout << std::endl;
        return false;*/

//c2v update
/*for (int ii = 0; ii < h_ins.check_num; ii++)
{
    cur_dc = h_ins.check_degreetable[ii];

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
        if (strcmp(check_op, "boxplus") == 0)
        {
            msg_c2v[h_ins.edge_c[ii][jj]] = check_node_operation(message);
        }
        else if (strcmp(check_op, "minsum") == 0)
        {
            msg_c2v[h_ins.edge_c[ii][jj]] = check_node_operation_min_double(message);
        }
        else if (strcmp(check_op, "att_minsum") == 0)
        {
            msg_c2v[h_ins.edge_c[ii][jj]] = check_node_operation_att_min_double(message, 0.7);
        }
        else if (strcmp(check_op, "minstar") == 0)
        {
            msg_c2v[h_ins.edge_c[ii][jj]] = check_node_operation_minstar(message);
        }
        else
        {
            std::cout << "Info: no such operation " << check_op << ".... please check again" << std::endl;
        }
    }
}

cur_dc = h_ins.check_degreetable[0];
for (int ii = 0; ii < cur_dc; ii++)
{
    std::cout << msg_c2v[h_ins.edge_c[0][ii]] << "  ";
}
std::cout << std::endl;
*/
