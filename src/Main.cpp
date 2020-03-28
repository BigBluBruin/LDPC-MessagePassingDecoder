
#include "ldpc_full_precision_decoder.h"


int main(int argc,char* argv[])
{
    ////---------------------------------------------------------------------------------------------------
    //argv[1]: check_op
    //argv[2]: eb_no
    //argv[3]: terget errors
    //argv[4]: max iteration time 
    //argv[5]: ind
    std::string parity_file="80211_irr_648_1296.txt";
    std::vector<double> parameters={std::stod(argv[2])};
    int target_error=std::stoi(argv[3]);
    ldpc_full_precision_decoder thisldpc(parity_file,parameters,target_error,std::stoi(argv[4]),argv[1]);
    thisldpc.fill_in();

    std::cout<<argc<<std::endl;
    if (argc==6)
    {
        thisldpc.main_simulation(argv[5]);
    }
    else if (argc==7)
    {
        thisldpc.main_simulation(argv[5],argv[6]);
    }

    //Let's track it
    // int total_num=35;
    // int iteration =0;
    // for (int cur_ind = 0; cur_ind<total_num; cur_ind++)
    // {
    //     std::string cur_initial_bit_filename="initial_bits_"+std::to_string(cur_ind)+".txt";
    //     std::string cur_cwd_filename = "wrong_cwd_"+std::to_string(cur_ind)+".txt";
    //     std::ifstream filehandle;
    //     std::vector<int> inital_bits(thisldpc.h_ins.vari_num,0);
    //     std::vector<int> final_bits(thisldpc.h_ins.vari_num,0);
    //     std::vector<double> cwds(thisldpc.h_ins.vari_num,0);
    //     filehandle.open(cur_initial_bit_filename);
    //     if(filehandle.is_open())
    //     {
    //         for(int ii=0; ii<thisldpc.h_ins.vari_num;ii++)
    //         {
    //             filehandle>>inital_bits[ii];
    //         }
    //     }
    //     else
    //     {
    //         std::cout<<"Info: Could not open INITIAL BIT file: "<<cur_initial_bit_filename<<". please check your input"<<std::endl;
    //         return 0;
    //     }
    //     filehandle.close();
    //     filehandle.open(cur_cwd_filename);
    //     if (filehandle.is_open())
    //     {
    //         for (int ii = 0; ii < thisldpc.h_ins.vari_num; ii++)
    //         {
    //             filehandle >> cwds[ii];
    //         }
    //     }
    //     else
    //     {
    //         std::cout << "Info: Could not open WRONG CWDs file: " << cur_cwd_filename << ". please check your input" << std::endl;
    //         return 0;
    //     }
    //     double cur_para = pow(10, (-0.1 * parameters[0]) / (2.0 * thisldpc.h_ins.rate));
    //     thisldpc.decoder_track(cwds,iteration,final_bits,inital_bits,cur_para,cur_ind);
    // }

}