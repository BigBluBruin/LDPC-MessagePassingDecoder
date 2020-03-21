
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
    
}