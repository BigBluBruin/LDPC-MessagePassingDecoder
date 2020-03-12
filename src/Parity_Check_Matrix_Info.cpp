#include "Parity_Check_Matrix_Info.h"

Parity_Check_Matrix_Info::Parity_Check_Matrix_Info(std::string Filename)
{
    filename=Filename;
    std::cout << "Input File: " << filename << "\n";
}

bool Parity_Check_Matrix_Info::Read_Parity_Check_Matrix()
{
    std::ifstream file_bar;
    file_bar.open(filename.c_str());
    int dummy;

    if (file_bar.is_open())
    {
        file_bar >> vari_num;
        file_bar >> check_num;
        file_bar >> edge_num;
        file_bar >> vari_max;
        file_bar >> check_max;
        //std::cout<<"vari: "<<vari_num<<"check: "<<check_num<<std::endl;
        rate = ((double)(vari_num - check_num)) / (double)(vari_num);
        vari_degreetable.resize(vari_num, 0);
        check_degreetable.resize(vari_num, 0);
        //read variable node degree table
        for (int i = 0; i < vari_num; i++)
        {
            file_bar >> dummy >> vari_degreetable[i];
        }
        //read check degree node table
        for (int i = 0; i < check_num; i++)
        {
            file_bar >> dummy >> check_degreetable[i];
            // std::cout<<check_degreetable[i]<<std::endl;
        }
        // read edge relationships
        edge_relation.resize(2, std::vector<int>(edge_num, 0));
        for (int i = 0; i < edge_num; i++)
        {
            //first vector --> variable
            //second vector --> check
            file_bar >> edge_relation[0][i] >> edge_relation[1][i];
        }
        file_bar.close();
        //find v->c dege relationships
        edge_v.resize(vari_num);
        for (int ii = 0; ii < vari_num; ii++)
        {
            int cur_dv = vari_degreetable[ii];
            int counter = 0;
            edge_v[ii].resize(cur_dv, 0);
            for (int jj = 0; jj < edge_num; jj++)
            {
                if (edge_relation[0][jj] == ii)
                {
                    edge_v[ii][counter] = jj;
                    counter++;
                }
                if (counter == cur_dv)
                {
                    break;
                }
            }
        }
        //find c-v edge relationships
        edge_c.resize(check_num);
        for (int ii = 0; ii < check_num; ii++)
        {
            int cur_dc = check_degreetable[ii];
            int counter = 0;
            edge_c[ii].resize(cur_dc, 0);
            for (int jj = 0; jj < edge_num; jj++)
            {
                if (edge_relation[1][jj] == ii)
                {
                    edge_c[ii][counter] = jj;
                    counter++;
                }
                if (counter == cur_dc)
                {
                    break;
                }
            }
        }

        return true;
    }
    else
    {
        std::cout << "can't open the file, check file name\n";
        return false;
    }
}

bool Parity_Check_Matrix_Info::Read_G_Matrix(std::string G_filename)
{
    std::ifstream myfile(G_filename);
    if (myfile.is_open())
    {
        int total_non_zero;
        myfile >> total_non_zero;
        G_info.assign(2, std::vector<unsigned>(total_non_zero, -1));
        for (int ii = 0; ii < total_non_zero; ii++)
        {
            myfile >> G_info[0][ii] >> G_info[1][ii];
        }
        std::cout << "G_matrix: " << G_filename << " has been successfully read" << std::endl;
        return true;
    }
    else
    {
        std::cout << "can't open G matrix file :" << G_filename << ". Please check again" << std::endl;
        return false;
    }
}

